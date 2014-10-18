module Transrate

  class ReadMetrics

    attr_reader :fragments_mapping
    attr_reader :p_good_mapping
    attr_reader :bad
    attr_reader :supported_bridges
    attr_reader :has_run
    attr_reader :read_length

    def initialize assembly
      @assembly = assembly
      @mapper = Snap.new
      self.initial_values

      load_executables
      @read_length = 100
    end

    def load_executables
      @bam_splitter = get_bin_path 'bam-split'
      @bam_reader = get_bin_path 'bam-read'
    end

    def get_bin_path bin
      which_bin = Cmd.new("which #{bin}")
      which_bin.run
      if !which_bin.status.success?
        raise IOError.new("ReadMetrics: could not find #{bin} in path")
      end
      which_bin.stdout.split("\n").first
    end

    def run left, right, insertsize:200, insertsd:50, threads:8
      # check all read files exist
      [left, right].each do |readfile|
        raise IOError.new "Read file is nil" if readfile.nil?
        readfile.split(",").each do |file|
          unless File.exist? file
            raise IOError.new "ReadMetrics: read file does not exist: #{file}"
          end
        end
      end

      # estimate max read length
      @read_length = get_read_length(left, right)

      # map reads
      @mapper.build_index(@assembly.file, threads)
      bamfile = @mapper.map_reads(@assembly.file, left, right,
                                  insertsize: insertsize,
                                  insertsd: insertsd,
                                  threads: threads)
      @fragments = @mapper.read_count

      # classify bam file into valid and invalid alignments
      sorted_bam = "#{File.basename(bamfile, '.bam')}.merged.sorted.bam"
      readsorted_bam = "#{File.basename(bamfile, '.bam')}.valid.sorted.bam"
      merged_bam = "#{File.basename(bamfile, '.bam')}.merged.bam"
      unless File.exist? sorted_bam
        unless File.exist? readsorted_bam
          valid_bam, invalid_bam = split_bam bamfile
          readsorted_bam = Samtools.readsort_bam valid_bam
          File.delete valid_bam
        end
      end

      # pass valid alignments to eXpress for assignment
      # always have to run the eXpress command to load the results
      assigned_bam = assign_and_quantify readsorted_bam
      File.delete readsorted_bam if File.exist? readsorted_bam

      # merge the assigned alignments back with the invalid ones
      unless File.exist? sorted_bam
        unless File.exist? merged_bam
          Samtools.merge_bam(invalid_bam, assigned_bam,
                             merged_bam, threads=threads)

          File.delete invalid_bam
          File.delete assigned_bam
        end
        sorted_bam = Samtools.sort_bam(merged_bam, [4, threads].min)
        File.delete merged_bam
      end

      # analyse the final mappings
      analyse_read_mappings(sorted_bam, insertsize, insertsd, true)

      @has_run = true
    end

    def read_stats
      {
        :fragments => @fragments,
        :fragments_mapped => @fragments_mapped,
        :p_fragments_mapped => @p_fragments_mapped,
        :good_mappings => @good,
        :p_good_mapping => @p_good_mapping,
        :bad_mappings => @bad,
        :potential_bridges => @potential_bridges,
        :bases_uncovered => @bases_uncovered,
        :p_bases_uncovered => @p_bases_uncovered,
        :contigs_uncovbase => @contigs_uncovbase,
        :p_contigs_uncovbase => @p_contigs_uncovbase,
        :contigs_uncovered => @contigs_uncovered,
        :p_contigs_uncovered => @p_contigs_uncovered,
        :contigs_lowcovered => @contigs_lowcovered,
        :p_contigs_lowcovered => @p_contigs_lowcovered,
        :contigs_segmented => @contigs_segmented,
        :p_contigs_segmented => @p_contigs_segmented,
        :contigs_good => @contigs_good,
        :p_contigs_good => @p_contigs_good
      }
    end

    def get_read_length(left, right)
      count=0
      file = File.open(left.split(",").first)
      name = file.readline.chomp
      seq = file.readline.chomp
      na = file.readline.chomp
      qual = file.readline.chomp
      read_length = 0
      while name and count < 5000 # get max read length from first 5000 reads
        read_length = [read_length, seq.length].max
        name = file.readline.chomp rescue nil
        seq = file.readline.chomp rescue nil
        na = file.readline.chomp rescue nil
        qual = file.readline.chomp rescue nil
        count+=1
      end
      read_length
    end

    def split_bam bamfile
      base = File.basename(bamfile, '.bam')
      valid = "#{base}.valid.bam"
      invalid = "#{base}.invalid.bam"
      if !File.exist? valid
        cmd = "#{@bam_splitter} #{bamfile}"
        splitter = Cmd.new cmd
        splitter.run
        if !splitter.status.success?
          logger.warn "Couldn't split bam file: #{bamfile}" +
                      "\n#{splitter.stdout}\n#{splitter.stderr}"
        end
      end
      if !File.exist? valid
        logger.warn "Splitting failed to create valid bam: #{valid}"
      end
      [valid, invalid]
    end

    def assign_and_quantify bamfile
      express = Express.new
      results = express.run(@assembly, bamfile)
      analyse_expression results.expression
      results.align_samp
    end

    def analyse_expression express_output
      express_output.each_pair do |name, expr|
        contig = @assembly[name]
        coverage = expr[:eff_count] * @read_length / expr[:eff_len]
        @contigs_uncovered += 1 if coverage < 1
        @contigs_lowcovered += 1 if coverage < 10
        contig.coverage = coverage
        contig.eff_len = expr[:eff_len]
        contig.eff_count = expr[:eff_count]
        contig.tpm = expr[:tpm]
      end
    end

    def analyse_read_mappings bamfile, insertsize, insertsd, bridge=true
      if File.exist?(bamfile) && File.size(bamfile) > 0
        csv_output = "#{File.basename(@assembly.file)}_bam_info.csv"
        csv_output = File.expand_path(csv_output)

        analyse_bam bamfile, csv_output
        # open output csv file
        @potential_bridges = 0

        CSV.foreach(csv_output, :headers => true,
                                :header_converters => :symbol,
                                :converters => :all) do |row|
          populate_contig_data row
        end
        @bad = @fragments_mapped - @good
      else
        logger.warn "couldn't find bamfile: #{bamfile}"
      end
      @assembly.assembly.each_pair do |name, contig|
        @contigs_good += 1 if contig.score >= 0.5
      end
      update_proportions
    end

    def update_proportions
      nbases = @assembly.n_bases.to_f
      ncontigs = @assembly.size.to_f

      @p_bases_uncovered = @bases_uncovered / nbases
      @p_contigs_uncovbase = @contigs_uncovbase / ncontigs
      @p_contigs_uncovered = @contigs_uncovered / ncontigs
      @p_contigs_lowcovered = @contigs_lowcovered / ncontigs
      @p_contigs_segmented = @contigs_segmented / ncontigs
      @p_contigs_good = @contigs_good / ncontigs

      @p_good_mapping = @good.to_f / @fragments.to_f
      @p_fragments_mapped = @fragments_mapped / @fragments.to_f
    end

    def analyse_bam bamfile, csv_output
      if !File.exist?(csv_output)
        cmd = "#{@bam_reader} #{bamfile} #{csv_output}"
        reader = Cmd.new cmd
        reader.run
        if !reader.status.success?
          logger.warn "couldn't get information from bam file: #{bamfile}"
        end
      end
    end

    def populate_contig_data row
      contig = @assembly[row[:name]]
      scale = 0.7
      contig.p_seq_true = (row[:p_seq_true] - scale) * (1.0 / (1 - scale))
      contig.uncovered_bases = row[:bases_uncovered]
      @bases_uncovered += contig.uncovered_bases
      if row[:fragments_mapped] and row[:fragments_mapped] > 0
        contig.p_good = row[:good]/row[:fragments_mapped].to_f
      end
      contig.p_not_segmented = row[:p_not_segmented]
      if contig.p_not_segmented < 0.5
        @contigs_segmented += 1
      end
      contig.in_bridges = row[:bridges]
      contig.p_unique = row[:p_unique]
      if row[:bridges] > 1
        @potential_bridges += 1
      end
      @fragments_mapped += row[:fragments_mapped]
      @good += row[:good]
      if row[:bases_uncovered] > 0
        @contigs_uncovbase += 1
      end
    end

    def initial_values
      @fragments = 0
      @fragments_mapped = 0
      @good = 0
      @bad = 0
      @bases_uncovered = 0
      @contigs_uncovbase = 0 # any base cov < 1
      @contigs_uncovered = 0 # mean cov < 1
      @contigs_lowcovered = 0 # mean cov < 10
      @contigs_segmented = 0 # p_not_segmented < 0.5
      @contigs_good = 0
    end

  end # ReadMetrics

end # Transrate
