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
      @mapper = Bowtie2.new
      self.initial_values

      which_bam = Cmd.new('which bam-read')
      which_bam.run
      if !which_bam.status.success?
        raise RuntimeError.new("could not find bam-read in path")
      end
      @bam_reader = which_bam.stdout.split("\n").first
      @read_length = 100
    end

    def run left, right, insertsize:200, insertsd:50, threads:8,
            sensitivity: "very-sensitive"
      [left, right].each do |readfile|
        raise IOError.new "Read file is nil" if readfile.nil?
        readfile.split(",").each do |file|
          unless File.exist? file
            raise IOError.new "ReadMetrics: read file does not exist: #{file}"
          end
        end
      end
      get_read_length(left, right)
      @mapper.build_index @assembly.file
      samfile = @mapper.map_reads(@assembly.file, left, right,
                                  insertsize: insertsize,
                                  insertsd: insertsd,
                                  threads: threads,
                                  sensitivity: sensitivity)
      @num_pairs = @mapper.read_count
      analyse_coverage(samfile)
      analyse_read_mappings(@sortedbam, insertsize, insertsd, true)

      @pr_good_mapping = @good.to_f / @num_pairs.to_f
      @percent_mapping = @total.to_f / @num_pairs.to_f * 100.0
      @pc_good_mapping = @pr_good_mapping * 100.0

      @has_run = true
    end

    # Generate per-base and contig read coverage statistics.
    # Note that contigs less than 200 bases long are ignored in this
    # analysis.
    def analyse_coverage(samfile)
      bam, @sortedbam, index = Samtools.sam_to_sorted_indexed_bam samfile
      bcf_file = Samtools.bam_to_bcf(@sortedbam, @assembly.file)
      #check bcf file is not empty
      line_count = 38 + @assembly.assembly.length
      if File.size(bcf_file) < 50*1024*1024
        line_checker = Cmd.new "wc -l #{bcf_file}"
        line_checker.run
        line_count = line_checker.stdout.split.first.to_i
      end
      # 37 is number of header lines in bcf file
      if line_count > 37 + @assembly.assembly.length
        load_bcf(bcf_file, @assembly.assembly.size, @read_length) # call to C
        # load contig data while len of contig > 0
        # if contig info returns contig length == -1 then stop because
        #  there are fewer contigs in the bcf file compared to the assembly
        #  because some contigs had no reads aligned to them (trinity)
        total_coverage = 0
        total_length = 0
        total_eff_length = 0
        total_eff_variance = 0
        total_mapq = 0
        i = 0
        len = 1
        while (i < @assembly.assembly.size)
          len = get_len(i)
          if len > 200
            contig_name = Bio::FastaDefline.new(get_contig_name(i)).entry_id
            contig = @assembly[contig_name]
            u = get_uncovered_bases(i)
            @n_uncovered_bases += u
            contig.uncovered_bases = u
            if u > 0
              @n_uncovered_base_contigs+=1
            end
            m = get_low_mapq_bases(i)
            contig.low_uniqueness_bases = m
            @n_low_uniqueness_bases += m
            contig.mean_coverage = get_total_coverage(i)/len.to_f
            if contig.mean_coverage < 1
              @n_uncovered_contigs+=1
            end
            if contig.mean_coverage < 10
              @n_lowcovered_contigs+=1
            end
            t = get_total_mapq(i)
            total_mapq += t
            contig.mean_mapq = t / len.to_f
            contig.variance = get_coverage_variance(i)

            # totals
            contig.effective_mean = get_mean_effective_coverage(i)
            contig.effective_variance = get_effective_coverage_variance(i)

            total_coverage += get_total_coverage(i)
            total_length += len
            total_eff_length += (len-200)
            total_eff_variance += get_effective_coverage_variance(i)*(len-200)
          end
          i += 1
        end
        @mean_coverage = total_coverage / total_length.to_f
        @coverage_variance = total_eff_variance / total_eff_length.to_f
        @mean_mapq = total_mapq / total_length.to_f
        @p_uncovered_bases = @n_uncovered_bases / total_length.to_f

        @p_uncovered_base_contigs = @n_uncovered_base_contigs /
                                    @assembly.assembly.size.to_f # any base cov < 1
        @p_uncovered_contigs = @n_uncovered_contigs /
                               @assembly.assembly.size.to_f # mean cov < 1
        @p_lowcovered_contigs = @n_lowcovered_contigs /
                                @assembly.assembly.size.to_f # mean cov < 10
        @p_low_uniqueness_bases = @n_low_uniqueness_bases / total_length.to_f
        free_contigs() #call to C to free the memory now we're finished with it
      else
        logger.warn "error creating bcf file. only has #{line_count} lines."
      end
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
        :mean_coverage => @mean_coverage,
        :coverage_variance => @coverage_variance,
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
      while name and count < 5000
        read_length = [read_length, seq.length].max
        name = file.readline.chomp rescue nil
        seq = file.readline.chomp rescue nil
        na = file.readline.chomp rescue nil
        qual = file.readline.chomp rescue nil
        count+=1
      end
      @read_length = read_length
    end

    def analyse_read_mappings bamfile, insertsize, insertsd, bridge=true
      if File.exist?(bamfile) && File.size(bamfile) > 0
        csv_output = "#{File.basename(@assembly.file)}_bam_info.csv"
        csv_output = File.expand_path(csv_output)
        if !File.exist?(csv_output)
          cmd = "#{@bam_reader} #{bamfile} #{csv_output}"
          reader = Cmd.new cmd
          reader.run
        end

        if !reader.status.success?
          logger.warn "couldn't get information from bam file"
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
          if row[:reads_mapped] and row[:reads_mapped]>0
            contig.p_good = row[:good]/row[:reads_mapped].to_f
          end
          @total_bases += row[:bases]
          contig.in_bridges = row[:bridges]
          if row[:bridges] > 1
            @supported_bridges += 1
          end
          @total += row[:both_mapped]
          @both_mapped += row[:both_mapped]
          @properly_paired += row[:properpair]
          @good += row[:good]

        end
        @bad = @num_pairs - @good
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
