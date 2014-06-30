module Transrate

  class ReadMetrics

    require 'bettersam'
    require 'bio-samtools'

    attr_reader :total
    attr_reader :bad
    attr_reader :supported_bridges
    attr_reader :pr_good_mapping
    attr_reader :percent_mapping
    attr_reader :prop_expressed
    attr_reader :has_run

    def initialize assembly
      @assembly = assembly
      @mapper = Bowtie2.new
      self.initial_values
    end

    def run left, right, insertsize=200, insertsd=50
      @mapper.build_index @assembly.file
      @num_pairs = `wc -l #{left}`.split(/\s+/)[0].to_i/4
      samfile = @mapper.map_reads(@assembly.file, left, right,
                                  insertsize: insertsize,
                                  insertsd: insertsd)
      self.analyse_read_mappings(samfile, insertsize, insertsd)
      self.analyse_coverage(samfile)
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
        :n_uncovered_bases => @n_uncovered_bases,
        :p_uncovered_bases => @p_uncovered_bases,
        :n_uncovered_base_contigs => @n_uncovered_base_contigs,
        :p_uncovered_base_contigs => @p_uncovered_base_contigs,
        :n_uncovered_contigs => @n_uncovered_contigs,
        :p_uncovered_contigs => @p_uncovered_contigs,
        :n_lowcovered_contigs => @n_lowcovered_contigs,
        :p_lowcovered_contigs => @p_lowcovered_contigs,:n_uncovered_contigs => @n_uncovered_contigs
      }
    end

    def analyse_read_mappings samfile, insertsize, insertsd, bridge=true
      @bridges = {} if bridge
      realistic_dist = self.realistic_distance(insertsize, insertsd)
      if File.exists?(samfile) && File.size(samfile) > 0
        ls = BetterSam.new
        rs = BetterSam.new
        sam = File.open(samfile)
        line = sam.readline
        while line=~/^@/
          line = sam.readline
        end
        while line
          ls.parse_line(line)
          if ls.mate_unmapped?
            self.check_read_single(ls)
            line = sam.readline rescue nil
          else
            line2 = sam.readline rescue nil
            if line2
              rs.parse_line(line2)
              self.check_read_pair(ls, rs, realistic_dist)
            end
            line = sam.readline rescue nil
          end
        end
      else
        raise "samfile #{samfile} not found"
      end
    end

    def initial_values
      @num_pairs = 0
      @total = 0
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
      @unexpressed_contigs = 0
      @expressed_contigs = 0
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
        @both_mapped += 1
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
          if count > 1
            f << [start, finish, count]
            @supported_bridges += 1
          end
        end
      end
    end

    def analyse_coverage bamfile
      bam = Bio::DB::Sam.new(:bam => bamfile, :fasta => @assembly.file)
      # get per-base coverage and calculate mean,
      # identify zero-coverage bases
      @n_uncovered_bases = 0
      @n_uncovered_base_contigs = 0 # any base cov < 1
      @n_uncovered_contigs = 0 # mean cov < 1
      @n_lowcovered_contigs = 0 # mean cov < 10
      @assembly.each do |contig|
        cov = bam.chromosome_coverage(contig.name,
                                      0,
                                      contig.length)
        contig.coverage = cov
        zerocov = 0
        total = 0
        cov.each do |e|
          total += e
          zerocov += 1 if e < 1
        end
        mean = total / cov.length.to_f
        @n_uncovered_bases += zerocov
        @n_uncovered_base_contigs += 1 if zerocov > 0
        @n_uncovered_contigs += 1 if mean < 1
        @n_lowcovered_contigs += 1 if mean < 10
      end
      @p_uncovered_bases = @n_uncovered_bases / @assembly.n_bases.to_f
      @p_uncovered_base_contigs = @n_uncovered_base_contigs /
                                  @assembly.size.to_f
      @p_uncovered_contigs = @n_uncovered_contigs / @assembly.size.to_f
      @p_lowcovered_contigs = @n_lowcovered_contigs / @assembly.size.to_f
    end

  end # ReadMetrics

end # Transrate
