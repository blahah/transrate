require 'set'
require 'crb-blast'

module Transrate
  
  class ComparativeMetrics

    attr_reader :rbh_per_contig
    attr_reader :rbh_per_reference
    attr_reader :reciprocal_hits
    attr_reader :has_run
    attr_reader :reference_coverage

    def initialize assembly, reference, threads
      @assembly = assembly
      @reference = reference
      @threads = threads
    end

    def run
      self.reciprocal_best_blast
      self.ortholog_hit_ratio
      self.chimeras
      @collapse_factor = collapse_factor @ra.target_results
      @reciprocal_hits = @ra.size
      @rbh_per_reference = @reciprocal_hits.to_f / @reference.size.to_f
      @reference_coverage = @ortholog_hit_ratio * @rbh_per_reference
      @rbh_per_contig = @reciprocal_hits.to_f / @assembly.assembly.size.to_f
      @has_run = true
    end

    def comp_stats
      {
        :reciprocal_hits => @reciprocal_hits,
        :rbh_per_contig => @rbh_per_contig,
        :rbh_per_reference => @rbh_per_reference,
        :reference_coverage => @reference_coverage,
        :ortholog_hit_ratio => @ortholog_hit_ratio,
        :collapse_factor => @collapse_factor
      }
    end

    def reciprocal_best_blast
      @ra = CRB_Blast.new @assembly.file, @reference.file
      @ra.run 1e-5, @threads
    end

    # coverage of contigs that have reciprocal hits
    # divided by
    # number of reciprocal targets
    def ortholog_hit_ratio
      return @ortholog_hit_ratio unless @ortholog_hit_ratio.nil?
      
      targets = Hash.new
      @ra.reciprocals.each_pair do |key, list|
        list.each do |hit|
          targets[hit.target] ||= [] # if key doesn't exist add it with a []
          targets[hit.target] << hit
        end
      end

      total_coverage=0
      total_length=0
      targets.each_pair do |key, list|

        blocks = []
        target_length = 0
        list.each do |hit|
          target_length = hit.tlen
          target_length *= 3 if @ra.target_is_prot
          start, stop = [hit.qstart, hit.qend].minmax
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
        # print "key = #{key}\t"
        # print "length_of_coverage = #{length_of_coverage}\t"
        # print "target_length = #{target_length}\n"
        total_coverage += length_of_coverage
        total_length += target_length
      end
      @ortholog_hit_ratio = total_coverage / total_length.to_f
    end

    def chimeras
      return @potential_chimera_ratio unless @potential_chimera_ratio.nil?
      potential_chimeras = 0

      @ra.reciprocals.each_pair do |key, list|
        blocks = []
        list.each do |hit|
          start, stop = [hit.qstart, hit.qend].minmax
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
          end
        end
        # if any blocks now overlap then extend one block and remove 
        # the other
        blocks.each_with_index do |block_a,a|
          blocks.each_with_index do |block_b,b|
            if a != b
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
              # elsif o == 5 # no overlap
                # do nothing
              # elsif o == 6 # no overlap
                # do nothing
              end
            end
          end
        end

        blocks.delete_if {|x| x[0]==-1 && x[1]==-1}
        if blocks.length > 1
          potential_chimeras += 1
        end
      end

      @potential_chimera_ratio = potential_chimeras / @ra.reciprocals.length.to_f
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
