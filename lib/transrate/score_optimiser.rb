module Transrate

  # This class is currently only used to calculate the basic transrate score.
  # In future it will be used to automatically optimised the score by
  # taking the optimal subset of contigs.
  class ScoreOptimiser

    def initialize assembly, read_metrics
      @assembly = assembly
      read_stats = read_metrics.read_stats
      @total = read_stats[:fragments]
      @good = read_stats[:good_mappings]
      raw_score
      @score = Score.new @contig_score, @good, @total
    end

    def raw_score
      scores = @assembly.assembly.values.map{ |c| c.score }
      @contig_score = geomean scores
      @contig_score * (@good / @total.to_f)
    end

    # Calculate the geometric mean of an array of numbers
    def geomean x
      sum = 0.0
      x.each{ |v| sum += Math.log(v) }
      sum /= x.size
      Math.exp sum
    end

  end # ScoreOptimiser

end
