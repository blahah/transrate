require 'set'
require 'crb-blast'

module Transrate

  class ComparativeMetrics

    attr_reader :rbh_per_contig
    attr_reader :rbh_per_reference
    attr_reader :reciprocal_hits
    attr_reader :has_run
    attr_reader :reference_coverage
    attr_reader :comp_stats
    attr_reader :n_chimeras, :p_chimeras

    def initialize assembly, reference, threads
      @assembly = assembly
      @reference = reference
      @threads = threads
      @comp_stats = Hash.new
    end

    def run
      @crbblast = reciprocal_best_blast
      @ortholog_hit_ratio = ortholog_hit_ratio @crbblast
      @collapse_factor = collapse_factor @crbblast.target_results
      @reciprocal_hits = @crbblast.size
      @rbh_per_reference = @reciprocal_hits.to_f / @reference.size.to_f
      @reference_coverage = @ortholog_hit_ratio * @rbh_per_reference
      @p_contigs_with_recip = @crbblast.reciprocals.size / @assembly.size.to_f
      @n_contigs_with_recip = @crbblast.reciprocals.size
      @p_refs_with_recip = @n_refs_with_recip / @reference.size.to_f
      chimeras @crbblast
      self.run_comp_stats
      @has_run = true
    end

    def run_comp_stats
      @comp_stats[:reciprocal_hits] = @reciprocal_hits
      @comp_stats[:rbh_per_contig] = @rbh_per_contig
      @comp_stats[:p_contigs_with_recip] = @p_contigs_with_recip
      @comp_stats[:n_contigs_with_recip] = @n_contigs_with_recip
      @comp_stats[:p_refs_with_recip] = @p_refs_with_recip
      @comp_stats[:n_refs_with_recip] = @n_refs_with_recip
      @comp_stats[:rbh_per_reference] = @rbh_per_reference
      @comp_stats[:reference_coverage] = @reference_coverage
      @comp_stats[:ortholog_hit_ratio] = @ortholog_hit_ratio
      @comp_stats[:collapse_factor] = @collapse_factor
      @comp_stats[:n_chimeras] = @n_chimeras
      @comp_stats[:p_chimeras] = @p_chimeras
    end

    def reciprocal_best_blast
      crbblast = CRB_Blast.new @assembly.file, @reference.file
      crbblast.run(1e-5, @threads, true)
      crbblast
    end

    # coverage of contigs that have reciprocal hits
    # divided by
    # number of reciprocal targets
    def ortholog_hit_ratio crbblast
      return @ortholog_hit_ratio unless @ortholog_hit_ratio.nil?

      targets = Hash.new
      crbblast.reciprocals.each_pair do |key, list|
        list.each do |hit|
          targets[hit.target] ||= [] # if key doesn't exist add it with a []
          targets[hit.target] << hit
        end
      end
      @n_refs_with_recip = targets.size
      total_coverage=0
      total_length=0
      targets.each_pair do |key, list|
        blocks = []
        target_length = 0
        list.each do |hit|
          target_length = hit.tlen
          if crbblast.target_is_prot
            target_length *= 3
            start, stop = [hit.tstart*3, hit.tend*3].minmax
          else
            start, stop = [hit.tstart, hit.tend].minmax
          end
          if blocks.empty?
            blocks << [start, stop]
          else
            found=false
            blocks.each do |block|
              # if query overlaps with any block extend that block
              o = overlap(block[0], block[1], start, stop)
              if o == 0 # perfect overlap
                found=true
              elsif o == 1 # partial overlap
                block[0] = start
                found=true
              elsif o == 2 # partial overlap
                block[1] = stop
                found=true
              elsif o == 3 # full overlap
                block[0] = start
                block[1] = stop
                found=true
              # elsif o == 4 # full overlap
                # nothing
              # elsif o == 5 || o == 6 # no overlap

              end
            end
            if !found
              blocks << [start, stop]
            end
            # if any blocks now overlap then extend one block and remove
            # the other
          end
        end
        blocks.each_with_index do |block_a,a|
          blocks.each_with_index do |block_b,b|
            if a!=b
              o = overlap(block_a[0], block_a[1], block_b[0], block_b[1])
              if o == 0 # perfect overlap
                block_b[0]=-1
                block_b[1]=-1
              elsif o == 1 # partial overlap
                block_a[0] = block_b[0]
                block_b[0] = -1
                block_b[1] = -1
              elsif o == 2 # partial overlap
                block_a[1] = block_b[1]
                block_b[0] = -1
                block_b[1] = -1
              elsif o == 3 # full overlap
                block_a[0] = block_b[0]
                block_a[1] = block_b[1]
                block_b[0] = -1
                block_b[1] = -1
              elsif o == 4 # full overlap
                block_b[0] = -1
                block_b[1] = -1
              # elsif o == 5 || o == 6# no overlap
                # do nothing
              # elsif  # no overlap
                # do nothing
              end
            end
          end # each_with_index b
        end # each_with_index a
        # sum blocks to find total coverage
        length_of_coverage=0
        blocks.each do |block|
          if block[0] and block[1]
            if block[0]>=0 and block[1]>=0
              length_of_coverage += block[1] - block[0] + 1
            end
          else
            puts "error: key = #{key}, #{blocks}"
          end
        end
        cov = [0.25, 0.5, 0.75, 0.85, 0.95]
        @cov ||= [0, 0, 0, 0, 0]
        p = length_of_coverage / target_length.to_f
        cov.each_with_index do |c, i|
          if p >= c
            @cov[i] +=1
          end
        end
        cov.each_with_index do |p, i|
          @comp_stats["cov#{(100*p).to_i}".to_sym] = @cov[i]
          @comp_stats["p_cov#{(100*p).to_i}".to_sym] = @cov[i]/@reference.size.to_f
        end
        total_coverage += length_of_coverage
        total_length += target_length
      end
      return ortholog_hit_ratio = total_coverage / total_length.to_f
    end

    def chimeras crbblast
      @n_chimeras = 0
      crbblast.reciprocals.each_pair do |key, list|
        p = 0
        list.each_with_index do |a, i|
          list.each_with_index do |b, j|
            if j>i
              if a.target == b.target
                astart, astop = [a.tstart, a.tend].minmax
                bstart, bstop = [b.tstart, b.tend].minmax

                oa = overlap_amount(astart, astop, bstart, bstop)
                if oa > 0.75
                  p += 1
                end
              else
                astart, astop = [a.qstart, a.qend].minmax
                bstart, bstop = [b.qstart, b.qend].minmax

                oa = overlap_amount(astart, astop, bstart, bstop)
                if oa < 0.25
                  p += 1
                end
              end
            end
          end
        end
        if p/list.size.to_f >= 0.5
          @n_chimeras += 1
        end
      end
      @p_chimeras = @n_chimeras / crbblast.reciprocals.length.to_f
    end

    def overlap(astart, astop, bstart, bstop)
      if astart == bstart and astop == bstop
        return 0
      elsif astart < bstart
        if astop > bstart
          if astop > bstop
            return 4
          else
            return 2
          end
        else
          return 5 # no overlap
        end
      else
        if bstop > astart
          if bstop > astop
            return 3
          else
            return 1
          end
        else
          return 6 # no overlap
        end
      end
    end

    def overlap_amount(astart, astop, bstart, bstop)
      if astart == bstart and astop == bstop
        return 1
      elsif astart < bstart
        if astop > bstart
          if astop > bstop
            return (bstop-bstart+1)/(astop-astart+1).to_f # 4
          else
            return (astop-bstart+1)/(bstop-astart+1).to_f # 2
          end
        else
          return 0 # 5 no overlap
        end
      else
        if bstop > astart
          if bstop > astop
            return (astop-astart+1)/(bstop-bstart+1).to_f # 3
          else
            return (bstop-astart+1)/(astop-bstart+1).to_f # 1
          end
        else
          return 0 # 6 no overlap
        end
      end
    end

    def collapse_factor hits=nil
      return @collapse_factor unless @collapse_factor.nil?
      targets = {}
      hits.each_pair do |query, list|
        list.each do |hit|
          target = hit.target
          unless targets.has_key? target
            targets[target] = Set.new
          end
          targets[target] << query
        end
      end
      sum = targets.values.reduce(0.0){ |summer, val| summer += val.size }
      sum / targets.size
    end

  end # ComparativeMetrics

end # Transrate
