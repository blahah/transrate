require 'set'
require 'inline'

module Transrate

  class ContigMetrics

    attr_reader :gc, :gc_skew, :at_skew, :cpg
    attr_reader :linguistic_complexity, :wootton_federhen, :shannon, :zlib_comp
    attr_reader :bases_n, :proportion_n
    attr_reader :has_run

    def initialize assembly
      @assembly = assembly
      self.initial_values
      @has_run = false
    end

    def initial_values
      @gc = -1
      @gc_skew = -1
      @at_skew = -1
      @cpg = -1
      @bases_n = -1
      @proportion_n = -1
      @linguistic_complexity = -1
      @wootton_federhen = -1   #
      @shannon = -1
      @zlib_comp = -1
    end

    def run
      calculate
      @has_run = true
    end

    def calculate
      total = 0
      a = 0
      c = 0
      g = 0
      t = 0
      @bases_n = 0
      cg = 0
      lc = 0
      set = Set.new
      k = 6
      @assembly.assembly.each do |entry|
        seq = entry.seq
        total += seq.length
        (0..seq.length-1).each do |i|
          a += 1 if "A" == seq[i].upcase
          c += 1 if "C" == seq[i].upcase
          g += 1 if "G" == seq[i].upcase
          t += 1 if "T" == seq[i].upcase
          @bases_n += 1 if "N" == seq[i].upcase
          if i > 0
            if seq[i-1].upcase == "C" and seq[i].upcase == "G"
              cg += 1
            end
          end
        end
        comp = get_linguistic_complexity(seq, 6)
        lc += comp
      end
      @gc = (g + c) / (a + c + g + t).to_f
      @gc_skew = (g - c) / (g + c).to_f
      @at_skew = (a - t) / (a + t).to_f
      @cpg = cg.to_f / (c * g) * total
      @linguistic_complexity = lc / @assembly.assembly.size.to_f
      @proportion_n = @bases_n / total.to_f
    end

    def results
      return if !@has_run
      return {'gc' => @gc,
              'gc_skew' => @gc_skew,
              'at_skew' => @at_skew,
              'cpg' => @cpg,
              'bases_n' => @bases_n,
              'proportion_n' => @proportion_n,
              'linguistic_complexity' => @linguistic_complexity
             }
    end

    def get_linguistic_complexity seq, k
      d = 4 ** k
      set = Set.new
      (0..seq.length-k).each do |i|
        set << seq.slice(i,k).upcase # slice(start, length)
      end # count how many kmers in seq
      set.size / d.to_f
    end

  end

end