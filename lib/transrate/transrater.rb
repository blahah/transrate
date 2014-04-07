module Transrate

  class Transrater

    attr_reader :assembly
    attr_reader :read_metrics
    attr_reader :comparative_metrics

    def initialize assembly, reference, left=nil, right=nil, insertsize=nil, insertsd=nil
      @assembly  = assembly.is_a?(Assembly)  ? assembly  : Assembly.new(assembly)
      @reference = reference.is_a?(Assembly) ? reference : Assembly.new(reference)
      @read_metrics = ReadMetrics.new @assembly
      @comparative_metrics = ComparativeMetrics.new(@assembly, @reference)
    end

    def run left, right, insertsize=nil, insertsd=nil
      @assembly.run
      @read_metrics.run(left, right)
      @comparative_metrics.run
    end

    def assembly_score
      pg = Metric.new('pg', @read_metrics.pr_good_mapping, 0.0)
      rc = Metric.new('rc', @comparative_metrics.reference_coverage, 0.0)
      #pce = Metric.new('pce', @read_metrics.prop_expressed, 0.0)
      puts "pg: #{pg.score}, rc: #{rc.score}" #, pce: #{pce.score}"
      @score = DimensionReduce.dimension_reduce([pg, rc])
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
