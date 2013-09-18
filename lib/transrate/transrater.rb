module Transrate

  class Transrater

    def initialize path
      @assembly = Assembly.new path
      @read_metrics = ReadMetrics.new @assembly
      @comparative_metrics = ComparativeMetrics.new @assembly
    end

    def run
      @assembly.run
      @read_metrics.run
      @comparative_metrics.run
    end

  end # Transrater

end # Transrate