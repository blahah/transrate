module Transrate

  class Transrater

    def initialize assembly, reference=nil, left=nil, right=nil, insertsize=nil, insertsd=nil
      @assembly  = assembly.is_a?(Assembly)  ? assembly  : Assembly.new(assembly)
      @reference = reference and reference.is_a?(Assembly) ? reference : Assembly.new(reference)
      @read_metrics = ReadMetrics.new @assembly
      @comparative_metrics = ComparativeMetrics.new(@assembly, @reference) if reference
    end

    def run left=nil, right=nil, insertsize=nil, insertsd=nil
      assembly_metrics
      if left && right
        read_metrics left, right
      end
      comparative_metrics
    end

    def assembly_score
      metrics = []
      assembly_metrics
      metrics << Metric.new("orf_percent", @assembly.orf_percent, 0.0)
      metrics << Metric.new("n50", @assembly.N50, 0.0)
      metrics.delete_if {|m| m.score==nil}
      @score = DimensionReduce.dimension_reduce(metrics)
    end

    def assembly_metrics
      @assembly.run unless @assembly.has_run
      @assembly
    end

    def read_metrics left=nil, right=nil
      @read_metrics.run(left, right) unless @read_metrics.has_run
      @read_metrics
    end

    def comparative_metrics
      @comparative_metrics.run unless @comparative_metrics.has_run
      @comparative_metrics
    end

    def all_metrics left, right, insertsize=nil, insertsd=nil
      self.run(left, right, insertsize, insertsd)
      all = @assembly.basic_stats
      all.merge!(@read_metrics.read_stats)
      all.merge!(@comparative_metrics.comp_stats)
      all[:score] = @score
      all
    end
    
  end # Transrater

end # Transrate
