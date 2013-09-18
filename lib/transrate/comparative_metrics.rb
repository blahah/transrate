require 'set'

module Transrate
  
  class ComparativeMetrics

    def initialize assembly, reference
      @assembly = assembly
      @reference = reference
      @usearch = Usearch.new
    end

    def run
      rbu = self.reciprocal_best_ublast
      @ortholog_hit_ratio = self.ortholog_hit_ratio rbu
      @collapse_factor = self.collapse_factor @ra.l2r_hits
      @reciprocal_hits = rbu.size
      {
        :reciprocal_hits => @reciprocal_hits,
        :ortholog_hit_ratio => @ortholog_hit_ratio,
        :collapse_factor => @collapse_factor
      }
    end

    def reciprocal_best_ublast
      @ra = ReciprocalAnnotation.new @assembly, @reference
      @ra.run
    end

    def ortholog_hit_ratio rbu
      rbu.reduce(0.0){ |sum, hit| sum += hit.last.tcov.to_f } / rbu.size
    end

    def collapse_factor hits
      targets = {}
      hits.each_pair do |query, hit|
        unless targets.has_key? query
          targets[query] = Set.new
        end
        targets[query] << hit.target
      end
      targets.values.reduce(0.0){ |sum, val| sum += val.size } / targets.size
    end

  end # ComparativeMetrics

end # Transrate