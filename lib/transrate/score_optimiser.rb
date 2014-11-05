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
      read_stats = read_metrics.read_stats
      total = read_stats.fragments
      good = read_stats.good
      raw_score
      @score = Score.new @contig_score, good, total
    end

    def raw_score
      @contig_score = geomean assembly.assembly.values.map do |contig|
        contig.score
      end
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
