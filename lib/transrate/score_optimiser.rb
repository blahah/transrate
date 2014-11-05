module Transrate

  class Score

    def initialize contig_score, good, total
      @n = 0
      @contig_score = contig_score
      @good = good
      @total = total
    end

    def remove_contig contig_good, contig_score

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

    def setup_contigs
      @contigs = @assembly.assembly.values.sort_by |contig|
        contig.score
      end
    end

    def optimal_score
      cutoff = 0.0
      good = @good
      score = raw_score
      last_contig_score = 0.0
      setup_contigs
      contig = next_contig
      n = 1
      while contig
        cutoff = contig.score
        good -= contig.good
        score = new_score(good, contig, @assembly.size - n, cutoff)

        n += 1
        contig = next_contig
      end
    end

    def next_contig
      @contigs.next
    rescue
      nil
    end

    def new_score(good, last_removed, n, old_contig_score, old_score)
      gm_old = old_contig_score ** (n + 1)
      gm_new = (gm_old / cutoff) ** (1.0 / last_removed)
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
