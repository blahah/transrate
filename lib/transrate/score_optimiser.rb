module Transrate

  # This class is used to calculate the basic transrate score and the
  # alternative weighted score.
  # It is also used to automatically optimise the score by
  # calculating a cutoff that maximises the number of reads that map
  # while also minimising the number of low scoring contigs.
  class ScoreOptimiser

    require 'csv'

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

    # Modification of the raw score, but where the contig quality
    # part of the score is the mean of each contig score is
    # multiplied by its relative expression.
    def weighted_score
      scores = @assembly.assembly.values.map{ |c| c.score * c.tpm }
      @contig_weighted = scores.inject{ |sum, s| sum + s }.to_f / scores.size
      @contig_weighted * (@good / @total.to_f)
    end

    def optimal_score(prefix='assembly')
      return [@optimal, @cutoff] unless @optimal.nil?
      product = 0
      good = 0
      @assembly.assembly.each do |key, contig|
        product += Math.log(contig.score)
        good += contig.good
      end
      count = @assembly.size
      cutoffscores = {}
      contigs_sorted = {}
      @assembly.assembly.sort_by { |k,v| v.score }.each do |a,b|
        contigs_sorted[a] = b
      end

      contigs_sorted.each do |key, contig|
        product -= Math.log(contig.score)
        good -= contig.good
        count -= 1
        score = Math.exp(product / count) * (good/@fragments.to_f)
        cutoffscores[contig.score] = score
      end
      @optimal = 0
      @cutoff = 0
      out = CSV.open("#{prefix}_score_optimisation.csv", 'w')
      out << %w[cutoff assembly_score]
      cutoffscores.each do |c, score|
        out << [c.round(5), score.round(5)]
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
