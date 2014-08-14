module Transrate

  class ReadMetrics

    require 'bettersam'

    attr_reader :total
    attr_reader :bad
    attr_reader :supported_bridges
    attr_reader :pr_good_mapping
    attr_reader :percent_mapping
    attr_reader :prop_expressed
    attr_reader :has_run
    attr_reader :total_bases

    def initialize assembly
      @assembly = assembly
      @mapper = Bowtie2.new
      self.initial_values
    end

    def run left, right, insertsize:200, insertsd:50, threads:8, singletons: nil
      [left, right].each do |readfile|
        raise IOError.new "Read file is nil" if readfile.nil?
        readfile.split(",").each do |file|
          unless File.exist? file
            raise IOError.new "ReadMetrics read file does not exist: #{file}"
          end
        end
      end
      @mapper.build_index @assembly.file
      samfile = @mapper.map_reads(@assembly.file, left, right,
                                  insertsize: insertsize,
                                  insertsd: insertsd,
                                  threads: threads)
      @num_pairs = @mapper.read_count
      # check_bridges
      analyse_read_mappings(samfile, insertsize, insertsd, true, singletons)
      analyse_coverage(samfile)
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

    def analyse_read_mappings samfile, insertsize, insertsd, bridge=true, singletons
      @bridges = {} if bridge
      realistic_dist = self.realistic_distance(insertsize, insertsd)
      if File.exists?(samfile) && File.size(samfile) > 0
        ls = BetterSam.new
        rs = BetterSam.new
        sam = File.open(samfile)
        out = File.open(singletons,'a') if singletons
        line = sam.readline
        while line and line=~/^@/
          line = sam.readline rescue nil
          out.puts(line) if singletons
        end
        while line
          ls.parse_line(line)
          if !ls.read_paired?
            out.puts(line) if singletons
          end
          lchrom = @assembly[ls.chrom]
          lchrom.edit_distance += ls.edit_distance
          lchrom.bases_mapped += ls.length
          @edit_distance += ls.edit_distance
          @total_bases += ls.length
          if ls.mate_unmapped?
            self.check_read_single(ls)
          else
            line2 = sam.readline rescue nil
            if line2
              rs.parse_line(line2)
              rchrom = (rs.chrom == ls.chrom) ? lchrom : @assembly[rs.chrom]
              rchrom.edit_distance += rs.edit_distance
              rchrom.bases_mapped += rs.length
              @edit_distance += rs.edit_distance
              @total_bases += rs.length
              self.check_read_pair(ls, rs, realistic_dist)
              out.puts(line) if ls.read_unmapped? if singletons
              out.puts(line2) if rs.read_unmapped? if singletons
            end
          end
          line = sam.readline rescue nil
        end
        check_bridges
        out.close if singletons
      else
        raise "samfile #{samfile} not found"
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

    def realistic_distance insertsize, insertsd
      insertsize + (3 * insertsd)
    end

    def check_read_single ls

    end

    def check_read_pair ls, rs, realistic_dist
      return unless ls.primary_aln?
      @total += 1
      if ls.both_mapped?
        # reads are paired
        @both_mapped += 1 if ls.primary_aln?
        if ls.read_properly_paired?
          # mapped in proper pair
          @properly_paired += 1
          self.check_orientation(ls, rs)
        else
          # not mapped in proper pair
          @improperly_paired += 1
          if ls.chrom == rs.chrom
            # both on same contig
            @same_contig += 1
            self.check_overlap_plausibility(ls, rs)
          else
            self.check_fragment_plausibility(ls, rs, realistic_dist)
          end
        end
      end
    end

    def check_orientation ls, rs
      if ls.pair_opposite_strands?
        # mates in proper orientation
        @proper_orientation += 1
        @good += 1
      else
        # mates in wrong orientation
        @improper_orientation += 1
        @bad += 1
      end
    end

    def check_overlap_plausibility ls, rs
      if Math.sqrt((ls.pos - rs.pos) ** 2) < ls.seq.length
        # overlap is realistic
        @realistic_overlap += 1
        self.check_orientation(ls, rs)
      else
        # overlap not realistic
        @unrealistic_overlap+= 1
        @bad += 1
      end
    end

    def check_fragment_plausibility ls, rs, realistic_dist
      # mates on different contigs
      # are the mapping positions within a realistic distance of
      # the ends of contigs?
      ldist = [ls.pos, ls.seq.length - ls.pos].min
      rdist = [rs.pos, rs.seq.length - rs.pos].min
      if ldist + rdist <= realistic_dist
        # increase the evidence for this bridge
        key = [ls.chrom, rs.chrom].sort.join("<>").to_sym
        if @bridges.has_key? key
          @bridges[key] += 1
        else
          @bridges[key] = 1
        end
        @realistic_fragment += 1
        @good += 1
      else
        @unrealistic_fragment += 1
        @bad += 1
      end
    end

    def check_bridges
      @supported_bridges = 0
      CSV.open('supported_bridges.csv', 'w') do |f|
        @bridges.each_pair do |b, count|
          start, finish = b.to_s.split('<>')
          @assembly[start].in_bridges += 1
          @assembly[finish].in_bridges += 1
          if count > 1
            f << [start, finish, count]
            @supported_bridges += 1
          end
        end
      end
    end

    # Generate per-base and contig read coverage statistics.
    # Note that contigs less than 200 bases long are ignored in this
    # analysis.
    def analyse_coverage samfile
      bamfile, sorted, index = Samtools.sam_to_sorted_indexed_bam samfile
      # get per-base coverage and calculate mean,
      # identify zero-coverage bases
      n_over_200, tot_length, tot_coverage, tot_mapq = 0, 0, 0, 0
      tot_variance, tot_eff_length = 0, 0
      @assembly.each_with_coverage(sorted, @assembly.file) do |contig,
                                                               coverage,
                                                               mapq|
        next if contig.length < 200
        n_over_200 += 1
        tot_length += contig.length
        tot_coverage += contig.load_coverage(coverage)
        tot_eff_length += contig.effective_length
        tot_mapq += contig.load_mapq(mapq)
        tot_variance += contig.effective_variance * contig.effective_length
        @n_uncovered_bases += contig.uncovered_bases
        @n_uncovered_base_contigs += 1 if contig.uncovered_bases > 0
        @n_uncovered_contigs += 1 if contig.mean_coverage < 1
        @n_lowcovered_contigs += 1 if contig.mean_coverage < 10
        @n_low_uniqueness_bases += contig.low_uniqueness_bases
      end
      @mean_coverage = (tot_coverage / tot_length.to_f).round(2)
      @mean_mapq = (tot_mapq / tot_length.to_f).round(2)
      @p_uncovered_bases = @n_uncovered_bases / tot_length.to_f
      @p_uncovered_base_contigs = @n_uncovered_base_contigs / n_over_200.to_f
      @p_uncovered_contigs = @n_uncovered_contigs / n_over_200.to_f
      @p_lowcovered_contigs = @n_lowcovered_contigs / n_over_200.to_f
      @p_low_uniqueness_bases = @n_low_uniqueness_bases / tot_length.to_f
      @coverage_variance = tot_variance / tot_eff_length
    end

  end # ReadMetrics

end # Transrate
