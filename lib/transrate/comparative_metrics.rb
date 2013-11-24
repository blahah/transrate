require 'set'

module Transrate
  
  class ComparativeMetrics

    attr_reader :rbh_per_contig
    attr_reader :reciprocal_hits
    attr_reader :has_run

    def initialize assembly, reference
      @assembly = assembly
      @reference = reference
      @usearch = Usearch.new
    end

    def run
      rbu = self.reciprocal_best_ublast
      @ortholog_hit_ratio = self.ortholog_hit_ratio rbu
      @collapse_factor = self.collapse_factor @ra.r2l_hits
      @reciprocal_hits = rbu.size
      @rbh_per_contig = @reciprocal_hits.to_f / @assembly.assembly.size.to_f
      @has_run = true
    end

    def comp_stats
      {
        :reciprocal_hits => @reciprocal_hits,
        :rbh_per_contig => @rbh_per_contig,
        :ortholog_hit_ratio => @ortholog_hit_ratio,
        :collapse_factor => @collapse_factor
      }
    end

    def reciprocal_best_ublast
      @ra = ReciprocalAnnotation.new @assembly, @reference
      @ra.run
    end

    def ortholog_hit_ratio rbu=nil
      return @ortholog_hit_ratio unless @ortholog_hit_ratio.nil? 
      rbu.reduce(0.0){ |sum, hit| sum += hit.last.tcov.to_f } / rbu.size
    end

    def collapse_factor hits=nil
      return @collapse_factor unless @collapse_factor.nil?
      targets = {}
      hits.each_pair do |query, hit|
        target = hit.target
        unless targets.has_key? target 
          targets[target] = Set.new
        end
        targets[target] << query
      end
      sum = targets.values.reduce(0.0){ |summer, val| summer += val.size }
      sum / targets.size
    end

  end # ComparativeMetrics

end # Transrate
