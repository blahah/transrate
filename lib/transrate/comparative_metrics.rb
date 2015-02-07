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

    def initialize assembly, reference, threads
      @assembly = assembly
      @reference = reference
      @threads = threads
      @comp_stats = Hash.new
    end

    def run
      @crbblast = reciprocal_best_blast
      @reference_coverage = coverage @crbblast
      @reciprocal_hits = @crbblast.size
      @rbh_per_reference = @reciprocal_hits.to_f / @reference.size.to_f
      @p_contigs_with_recip = @crbblast.reciprocals.size / @assembly.size.to_f
      @n_contigs_with_recip = @crbblast.reciprocals.size
      count_ref_crbbs
      @p_refs_with_recip = @n_refs_with_recip / @reference.size.to_f
      self.run_comp_stats
      @has_run = true
    end

    def run_comp_stats
      @comp_stats[:CRBB_hits] = @reciprocal_hits # CRBB hits
      @comp_stats[:p_contigs_with_CRBB] = @p_contigs_with_recip
      @comp_stats[:n_contigs_with_CRBB] = @n_contigs_with_recip
      @comp_stats[:p_refs_with_CRBB] = @p_refs_with_recip
      @comp_stats[:n_refs_with_CRBB] = @n_refs_with_recip
      @comp_stats[:rbh_per_reference] = @rbh_per_reference
      @comp_stats[:reference_coverage] = @reference_coverage
    end

    def reciprocal_best_blast
      crbblast = CRB_Blast::CRB_Blast.new @assembly.file, @reference.file
      crbblast.run(1e-5, @threads, true)
      crbblast
    end

    # coverage of contigs that have reciprocal hits
    # divided by number of reciprocal targets
    def coverage crbblast
      return @reference_coverage unless @reference_coverage.nil?
      crbblast.reciprocals.each do |key, list|
        list.each_with_index do |hit, i|
          unless @reference.assembly.key? hit.target
            raise TransrateError.new "#{hit.target} not in reference"
          end
          @reference[hit.target].hits << hit

          unless @assembly.assembly.key? hit.query
            raise TransrateError.new "#{hit.query} not in assembly"
          end
          contig = @assembly[hit.query]
          contig.has_crb = true
          # how much of the reference is covered by this single contig
          if crbblast.target_is_prot
            contig.reference_coverage =
                        (hit.alnlen - hit.mismatches - hit.gaps) / (3*hit.tlen)
          else
            contig.reference_coverage =
                            (hit.alnlen - hit.mismatches - hit.gaps) / hit.tlen
          end
          contig.hits << hit
        end
      end
      total_coverage = 0
      total_length = 0
      cov = [0.25, 0.5, 0.75, 0.85, 0.95]
      @cov ||= [0, 0, 0, 0, 0]
      @reference.each_value do |ref_contig|
        key = ref_contig.name
        list = ref_contig.hits
        if crbblast.target_is_prot
          total_length += ref_contig.length * 3
        else
          total_length += ref_contig.length
        end
        next if list.empty?
        blocks = []
        target_length = 0
        list.each do |hit|
          target_length = hit.tlen
          if crbblast.target_is_prot
            target_length *= 3
            start, stop = [hit.tstart, hit.tend].minmax
            start = start*3-2
            stop = stop*3
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
              elsif o == 4 # full overlap
                found=true
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
        length_of_coverage = calculate_coverage blocks
        if target_length > 0
          ref_p = length_of_coverage / target_length.to_f
        else
          ref_p = 0
        end
        ref_contig.reference_coverage = ref_p

        cov.each_with_index do |c, i|
          if ref_p >= c
            @cov[i] +=1
          end
        end

        total_coverage += length_of_coverage
      end

      cov.each_with_index do |p, i|
        @comp_stats["cov#{(100*p).to_i}".to_sym] = @cov[i]
        @comp_stats["p_cov#{(100*p).to_i}".to_sym] =
          @cov[i]/@reference.size.to_f
      end
      total_coverage / total_length.to_f
    end

    # Calculate the total coverage from a set of coverage blocks
    def calculate_coverage blocks
      coverage = 0
      blocks.each do |block|
        if block[0] and block[1]
          if block[0]>=0 and block[1]>=0
            coverage += block[1] - block[0] + 1
          end
        else
          puts "error: key = #{key}, #{blocks}"
        end
      end
      coverage
    end

    # Count reference proteins with at least one recprocal hit
    def count_ref_crbbs
      @n_refs_with_recip = @reference.assembly.inject(0) do |sum, entry|
        name, contig = entry
        sum + (contig.hits.length > 0 ? 1 : 0)
      end
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

  end # ComparativeMetrics

end # Transrate
