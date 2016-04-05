require 'set'

module Transrate

  class ContigMetrics

    attr_reader :gc_prop
    attr_reader :bases_n, :proportion_n
    attr_reader :has_run

    def initialize assembly
      @assembly = assembly
      self.initial_values
      @has_run = false
    end

    def initial_values
      @gc_prop = -1
      @bases_n = -1
      @proportion_n = -1
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
      @assembly.assembly.each_value do |contig|
        total += contig.length
        a += contig.bases_a
        c += contig.bases_c
        g += contig.bases_g
        t += contig.bases_t
        @bases_n += contig.bases_n
      end
      @gc_prop = (g + c) / (a + c + g + t).to_f
      @proportion_n = @bases_n / total.to_f
    end

    def results
      return if !@has_run
      return {
        'gc' => @gc_prop,
        'bases_n' => @bases_n,
        'proportion_n' => @proportion_n,
       }
    end

  end

end
