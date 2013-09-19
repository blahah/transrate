module Transrate

  class Transrater

    def initialize assembly, reference, left, right, insertsize=nil, insertsd=nil
      @assembly = Assembly.new assembly
      @reference = Assembly.new reference
      @read_metrics = ReadMetrics.new @assembly
      @comparative_metrics = ComparativeMetrics.new(@assembly, @reference)
      self.run(left, right, insertsize, insertsd)
    end

    def run left, right, insertsize=nil, insertsd=nil
      @assembly.run
      @read_metrics.run(left, right)
      @comparative_metrics.run
    end

  end # Transrater

end # Transrate