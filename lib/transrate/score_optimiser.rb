module Transrate

  # This class is currently only used to calculate the basic transrate score.
  # In future it will be used to automatically optimised the score by
  # taking the optimal subset of contigs.
  class ScoreOptimiser

    def initialize assembly, read_metrics
      @assembly = assembly
      @fragments = read_metrics.fragments
      read_stats = read_metrics.read_stats
      @total = read_stats[:fragments]
      @good = read_stats[:good_mappings]
    end

    def raw_score
      scores = @assembly.assembly.values.map{ |c| c.score }
      @contig_score = geomean scores
      @contig_score * (@good / @total.to_f)
    end

    def optimal_score
      return [@optimal, @cutoff] unless @optimal.nil?
      product = 0
      good = 0
      @assembly.assembly.each do |key, contig|
        product += Math.log(contig.score)
        good += contig.good
      end
      count = @assembly.size
      cutoffscores = {}
      contigs_sorted = @assembly.assembly.sort_by {|k,v| v.score}.to_h

      contigs_sorted.each do |key, contig|
        product -= Math.log(contig.score)
        good -= contig.good
        count -= 1
        score = Math.exp(product / count) * (good/@fragments.to_f)
        cutoffscores[contig.score] = score
      end
      @optimal = 0
      @cutoff = 0
      cutoffscores.each do |c, score|
        if score > @optimal
          @optimal = score
          @cutoff = c
        end
      end
      return [@optimal, @cutoff]
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
