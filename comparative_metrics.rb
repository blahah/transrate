require_relative 'reciprocal_annotation'
require 'set'

class ComparativeMetrics

  def initialize assembly, reference
    @assembly = assembly
    @reference = reference
    @usearch = Usearch.new
  end

  def run
    rbu = self.reciprocal_best_ublast
    ohr = self.ortholog_hit_ratio rbu
    cf = self.collapse_factor rbu.l2r_hits
    {
      :rbu => rbu,
      :ohr => ohr,
      :cf => cf
    }
  end

  def reciprocal_best_ublast
    ra = ReciprocalAnnotation.new @assembly, @reference
    ra.run
  end

  def ortholog_hit_ratio rbu
    rbu.reduce(0.0){ |sum, hit| sum += hit.tcov } / rbu.size
  end

  def collapse_factor hits
    targets = Hash.new(Set.new)
    hits.each_pair do |query, hit|
      targets[query] << hit.target
    end
    targets.each_pair.reduce(0.0){ |sum, key, val| sum += val.size } / targets.size
  end

end # ComparativeMetrics