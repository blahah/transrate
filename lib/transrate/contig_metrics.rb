require 'set'

module Transrate

  class ContigMetrics

    attr_reader :gc_prop, :gc_skew, :at_skew, :cpg_ratio
    attr_reader :linguistic_complexity, :wootton_federhen, :shannon, :zlib_comp
    attr_reader :bases_n, :proportion_n
    attr_reader :has_run

    def initialize assembly
      @assembly = assembly
      self.initial_values
      @has_run = false
    end

    def initial_values
      @gc_prop = -1
      @gc_skew = -1
      @at_skew = -1
      @cpg = -1
      @bases_n = -1
      @proportion_n = -1
      @linguistic_complexity = -1
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
      cpg_count = 0
      lc = 0
      k = 6
      @assembly.assembly.each_value do |contig|
        total += contig.length
        a += contig.bases_a
        c += contig.bases_c
        g += contig.bases_g
        t += contig.bases_t
        @bases_n += contig.bases_n
        cpg_count += contig.cpg_count
        lc += contig.linguistic_complexity k
      end
      @gc_prop = (g + c) / (a + c + g + t).to_f
      @linguistic_complexity = lc / @assembly.assembly.size.to_f
      @proportion_n = @bases_n / total.to_f
    end

    def results
      return if !@has_run
      return {'gc' => @gc_prop,
              'gc_skew' => @gc_skew,
              'at_skew' => @at_skew,
              'cpg_ratio' => @cpg_ratio,
              'bases_n' => @bases_n,
              'proportion_n' => @proportion_n,
              'linguistic_complexity' => @linguistic_complexity
             }
    end

  end

end
