module Transrate

  class DimensionReduce

    def self.dimension_reduce(metrics)
      total = 0
      metrics.each do |metric|
        o = metric.origin
        w = metric.weighting
        a = metric.score
        total += w * ((o - a) ** 2)
      end
      Math.sqrt(total) / metrics.length
    end
      
  end # DimensionReduce

end # Transrate
