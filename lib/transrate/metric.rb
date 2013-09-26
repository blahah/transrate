module Transrate

  class Metric

    attr_reader :origin, :score, :name, :weighting
    
    def initialize(name, score, origin)
      @origin = origin
      @score = score
      @name = name
      @weighting = 1
    end

  end # Metric

end # Transrate
