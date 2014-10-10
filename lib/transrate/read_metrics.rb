module Transrate

  class ReadMetrics

    attr_reader :total
    attr_reader :bad
    attr_reader :supported_bridges
    attr_reader :pr_good_mapping
    attr_reader :percent_mapping
    attr_reader :prop_expressed
    attr_reader :has_run
    attr_reader :total_bases
    attr_reader :read_length

    def initialize assembly
      @assembly = assembly
      @mapper = Snap.new
      self.initial_values

      which_bam = Cmd.new('which bam-read')
      which_bam.run
      if !which_bam.status.success?
        raise RuntimeError.new("could not find bam-read in path")
      end
      @bam_reader = which_bam.stdout.split("\n").first
      @read_length = 100
      @mean_mapq = 0
    end

    def run left, right, insertsize:200, insertsd:50, threads:8
      [left, right].each do |readfile|
        raise IOError.new "Read file is nil" if readfile.nil?
        readfile.split(",").each do |file|
          unless File.exist? file
            raise IOError.new "ReadMetrics: read file does not exist: #{file}"
          end
        end
      end
      get_read_length(left, right)
      @mapper.build_index(@assembly.file, threads)
      bamfile = @mapper.map_reads(@assembly.file, left, right,
                                  insertsize: insertsize,
                                  insertsd: insertsd,
                                  threads: threads)
      @num_pairs = @mapper.read_count
      # analyse_coverage(samfile)
      analyse_read_mappings(bamfile, insertsize, insertsd, true)

      @pr_good_mapping = @good.to_f / @num_pairs.to_f
      @percent_mapping = @total.to_f / @num_pairs.to_f * 100.0
      @pc_good_mapping = @pr_good_mapping * 100.0

      @has_run = true
    end

    def read_stats
      {
        :num_pairs => @num_pairs,
        :total_mappings => @total,
        :percent_mapping => @percent_mapping,
        :good_mappings => @good,
        :pc_good_mapping => @pc_good_mapping,
        :bad_mappings => @bad,
        :potential_bridges => @supported_bridges,
        :mean_mapq => @mean_mapq,
        :n_uncovered_bases => @n_uncovered_bases,
        :p_uncovered_bases => @p_uncovered_bases,
        :n_uncovered_base_contigs => @n_uncovered_base_contigs,
        :p_uncovered_base_contigs => @p_uncovered_base_contigs,
        :n_uncovered_contigs => @n_uncovered_contigs,
        :p_uncovered_contigs => @p_uncovered_contigs,
        :n_lowcovered_contigs => @n_lowcovered_contigs,
        :p_lowcovered_contigs => @p_lowcovered_contigs,
        :edit_distance_per_base => @edit_distance / @total_bases.to_f,
        :n_low_uniqueness_bases => @n_low_uniqueness_bases,
        :p_low_uniqueness_bases => @p_low_uniqueness_bases
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
      @read_length = read_length
    end

    def analyse_expression express_output
      # set n_uncovered_contigs here. unexpressed contigs
      # @n_uncovered_contigs += 1 if expression < 1
    end

    def analyse_read_mappings bamfile, insertsize, insertsd, bridge=true
      if File.exist?(bamfile) && File.size(bamfile) > 0
        csv_output = "#{File.basename(@assembly.file)}_bam_info.csv"
        csv_output = File.expand_path(csv_output)
        if !File.exist?(csv_output)
          cmd = "#{@bam_reader} #{bamfile} #{csv_output}"
          reader = Cmd.new cmd
          reader.run
          if !reader.status.success?
            logger.warn "couldn't get information from bam file"
          end
        end
        # open output csv file
        @supported_bridges = 0

        CSV.foreach(csv_output, :headers => true,
                                :header_converters => :symbol,
                                :converters => :all) do |row|
          contig = @assembly[row[:name]]
          @edit_distance += row[:edit_distance]
          contig.edit_distance = row[:edit_distance]
          contig.bases_mapped = row[:bases]
          contig.uncovered_bases = row[:bases_uncovered]

          if row[:reads_mapped] and row[:reads_mapped]>0
            contig.p_good = row[:good]/row[:reads_mapped].to_f
          end
          contig.p_not_segmented = row[:p_not_segmented]
          @total_bases += row[:bases]
          contig.in_bridges = row[:bridges]
          if row[:bridges] > 1
            @supported_bridges += 1
          end
          @total += row[:both_mapped]
          @both_mapped += row[:both_mapped]
          @properly_paired += row[:properpair]
          @good += row[:good]
          if row[:bases_uncovered]>0
            @n_uncovered_base_contigs += 1
          end
          @mean_mapq += row[:mapq]
        end
        @bad = @num_pairs - @good
        @mean_mapq /= @assembly.length
      else
        logger.warn "couldn't find bamfile"
      end

    end

    def initial_values
      @num_pairs = 0
      @total = 0
      @total_bases = 0
      @good = 0
      @bad = 0
      @both_mapped = 0
      @properly_paired = 0
      @improperly_paired = 0
      @proper_orientation = 0
      @improper_orientation = 0
      @same_contig = 0
      @realistic_overlap = 0
      @unrealistic_overlap = 0
      @realistic_fragment = 0
      @unrealistic_fragment = 0
      @n_uncovered_bases = 0
      @n_uncovered_base_contigs = 0 # any base cov < 1
      @n_uncovered_contigs = 0 # mean cov < 1
      @n_lowcovered_contigs = 0 # mean cov < 10
      @edit_distance = 0
      @n_low_uniqueness_bases = 0
    end

  end # ReadMetrics

end # Transrate
