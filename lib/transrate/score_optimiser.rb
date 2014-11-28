module Transrate

  class Score

    def initialize contig_score, good, total
      @n = 0
      @contig_score = contig_score
      @good = good
      @total = total
    end

  end # Score

  class ScoreOptimiser

    def initialize assembly, read_metrics
      @assembly = assembly
      total = read_metrics.fragments
      good = read_metrics.good
      raw_score
      @score = Score.new @contig_score, good, total
    end

    def raw_score
      contig_scores = @assembly.assembly.values.map do |contig|
        contig.score
      end
      @contig_score = geomean contig_scores
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
