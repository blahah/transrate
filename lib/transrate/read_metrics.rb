module Transrate

  class ReadMetricsError < TransrateError; end

  class ReadMetrics

    attr_reader :fragments, :fragments_mapping, :p_good_mapping
    attr_reader :good, :bad
    attr_reader :supported_bridges
    attr_reader :has_run
    attr_reader :read_length

    def initialize assembly
      @assembly = assembly # Transrate::Assembly
      @mapper = Snap.new
      @salmon = Salmon.new
      self.initial_values

      load_executables
      @read_length = 100
    end

    def load_executables
      @bam_reader = get_bin_path 'bam-read'
    end

    def get_bin_path bin
      which_bin = Cmd.new("which #{bin}")
      which_bin.run
      if !which_bin.status.success?
        raise TransrateIOError.new("ReadMetrics: could not find #{bin} in path")
      end
      which_bin.stdout.split("\n").first
    end

    def run left, right, threads:8
      # check all read files exist
      [left, right].each do |readfile|
        raise TransrateIOError.new "Read file is nil" if readfile.nil?
        readfile.split(",").each do |file|
          unless File.exist? file
            raise TransrateIOError.new "ReadMetrics: read file does not " +
                                       "exist: #{file}"
          end
        end
      end

      # estimate max read length
      @read_length = get_read_length(left, right)

      # map reads
      @mapper.build_index(@assembly.file, threads)
      bamfile = @mapper.map_reads(@assembly.file, left, right,
                                  threads: threads)
      @fragments = @mapper.read_count

      assigned_bam = "postSample.bam"
      final_bam = "#{File.basename(bamfile, '.bam')}.assigned.bam"

      # check for latest files first and create what is needed
      if !File.exist?(final_bam)
        if !File.exist?(assigned_bam)
          assigned_bam = assign_and_quantify(bamfile, threads)
        end
        if File.exist?(assigned_bam)
          File.rename(assigned_bam, final_bam)
        else
          logger.error "Couldn't find #{assigned_bam} to rename"
          raise ReadMetricsError
        end
      end
      # analyse the final mappings
      analyse_read_mappings final_bam

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
        :p_contigs_segmented => @p_contigs_segmented
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

    def assign_and_quantify(bamfile, threads)
      @salmon.run(@assembly, bamfile, threads)
    end

    def analyse_expression salmon_output
      salmon_output.each_pair do |name, expr|
        contig_name = Bio::FastaDefline.new(name.to_s).entry_id
        contig_name.gsub!(/;$/, '') # trim trailing semicolon
        contig = @assembly[contig_name]
        if expr[:eff_len]==0
          coverage = 0
        else
          coverage = expr[:eff_count] * @read_length / expr[:eff_len]
        end
        @contigs_uncovered += 1 if coverage < 1
        @contigs_lowcovered += 1 if coverage < 10
        contig.coverage = coverage.round(2)
        contig.eff_length = expr[:eff_len]
        contig.eff_count = expr[:eff_count]
        contig.tpm = expr[:tpm]
      end
    end

    def analyse_read_mappings bamfile
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
        raise TransrateError.new "couldn't find bamfile: #{bamfile}"
      end
      salmon_results = "#{File.basename @assembly.file}_quant.sf"

      if File.exist?(salmon_results)
        analyse_expression(@salmon.load_expression(salmon_results))
      else
        abort "Can't find #{salmon_results}"
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

      @p_good_mapping = @good.to_f / @fragments.to_f
      @p_fragments_mapped = @fragments_mapped / @fragments.to_f
    end

    def analyse_bam bamfile, csv_output
      if !File.exist?(csv_output)
        cmd = "#{@bam_reader} #{bamfile} #{csv_output} 0.7"
        reader = Cmd.new cmd
        reader.run
        if !reader.status.success?
          msg = "Couldn't get information from bam file: #{bamfile}\n"
          msg << "#{reader.stdout}\n#{reader.stderr}"
          raise TransrateError.new msg
        end
      end
    end

    def populate_contig_data row
      name = Bio::FastaDefline.new(row[:name].to_s).entry_id
      name.gsub!(/;$/, '') # trim trailing semicolon
      contig = @assembly[name]
      contig.p_seq_true = row[:p_seq_true]
      contig.uncovered_bases = row[:bases_uncovered]
      @bases_uncovered += contig.uncovered_bases
      if row[:fragments_mapped] and row[:fragments_mapped] > 1
        contig.p_good = row[:good]/row[:fragments_mapped].to_f
      end
      contig.p_not_segmented = row[:p_not_segmented]
      if contig.p_not_segmented < 0.5
        @contigs_segmented += 1
      end
      contig.in_bridges = row[:bridges]
      if row[:bridges] > 1
        @potential_bridges += 1
      end
      @fragments_mapped += row[:fragments_mapped]
      contig.good = row[:good]
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
    end

  end # ReadMetrics

end # Transrate
