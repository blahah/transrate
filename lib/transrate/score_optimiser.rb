module Transrate

  # This class is currently only used to calculate the basic transrate score.
  # In future it will be used to automatically optimised the score by
  # taking the optimal subset of contigs.
  class ScoreOptimiser

    def initialize assembly, read_metrics
      @assembly = assembly
      read_stats = read_metrics.read_stats
      total = read_stats[:fragments]
      good = read_stats[:good_mappings]
      raw_score
      @score = Score.new @contig_score, good, total
    end

    def raw_score
      @contig_score = geomean @assembly.assembly.values.map do |contig|
        contig.score
      end
      @contig_score * (@good / @total.to_f)
    end

  end # ScoreOptimiser

end
