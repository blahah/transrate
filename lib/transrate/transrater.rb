module Transrate

  class Transrater

    attr_reader :assembly
    attr_reader :read_metrics
    attr_reader :comparative_metrics

    def initialize assembly, reference, left, right, insertsize=nil, insertsd=nil
      @assembly = assembly.is_a?(Assembly) ? assembly : Assembly.new(assembly)
      @reference = reference.is_a?(Assembly) ? reference : Assembly.new(reference)
      @read_metrics = ReadMetrics.new @assembly
      @comparative_metrics = ComparativeMetrics.new(@assembly, @reference)
      self.run(left, right, insertsize, insertsd)
    end

    def run left, right, insertsize=nil, insertsd=nil
      @assembly.run
      @read_metrics.run(left, right)
      @comparative_metrics.run
    end

    def assembly_score
      # percent of reads mapping in biologically feasible way
      pg = @read_metrics.pc_good_mapping
      pg = Metric.new('pg', pg, 0)
      # number of reciprocal hits
      rh = @comparative_metrics.reciprocal_hits
      rh = Metric.new('rh', rh, 0)
      # ortholog hit ratio
      ohr = @comparative_metrics.ortholog_hit_ratio
      ohr = Metric.new('ohr', ohr, 0)
      # TODO: add number of expressed transcripts
      # TODO: consider using proportion of contigs
      #       that have reciprocal hits?
      DimensionReduce.dimension_reduce([pg, rh])
    end
    
  end # Transrater

end # Transrate
