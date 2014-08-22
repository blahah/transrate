module Transrate

  # A transrater runs all types of metrics on an assembly.
  #
  # @!attribute [r] assembly
  #   @return [Assembly, String] an Assembly or the path to an assembly
  # @!attribute [r] read_metrics
  #   @return [Hash] the read metrics if they have been calculated
  # @!attribute [r] comparative_metrics
  #   @return [hash] the comparative metrics if they have been calculated
  class Transrater

    attr_reader :assembly
    attr_reader :read_metrics

    # A new Transrater
    #
    # @param assembly [Assembly, String] the Assembly or path to the FASTA
    # @param reference [Assembly, String] the reference Assembly or
    #   path to the FASTA
    # @param left [String] path to the left reads
    # @param right [String] path to the right reads
    # @param insertsize [Integer] mean insert size of the read pairs
    # @param insertsd [Integer] standard deviation of the read pair insert size
    def initialize(assembly, reference,
                   left: nil, right: nil,
                   insertsize: nil, insertsd: nil,
                   threads: 1)
      if assembly
        if assembly.is_a?(Assembly)
          @assembly = assembly
        else
          @assembly = Assembly.new(assembly)
        end
        @read_metrics = ReadMetrics.new @assembly
      else
        raise RuntimeError.new("assembly is nil")
      end

      if reference
        if reference.is_a?(Assembly)
          @reference = reference
        else
          @reference = Assembly.new(reference)
        end
        @comparative_metrics = ComparativeMetrics.new(@assembly,
                                                      @reference,
                                                      threads)
      end
      @threads = threads
    end

    # Run all analyses
    #
    # @param left [String] path to the left reads
    # @param right [String] path to the right reads
    # @param insertsize [Integer] mean insert size of the read pairs
    # @param insertsd [Integer] standard deviation of the read pair insert size
    def run left=nil, right=nil, insertsize=nil, insertsd=nil
      assembly_metrics
      if left && right
        read_metrics left, right
      end
      comparative_metrics
    end

    # Calculate the geometric mean of an array of numbers
    def geomean(x)
      sum = 0.0
      x.each{ |v| sum += Math.log(v) }
      sum /= x.size
      Math.exp(sum)
    end

    # Reduce all metrics for the assembly to a single quality score
    # by taking the geometric mean of the scores for all contigs
    #
    #
    # @return [Integer] the assembly score
    def assembly_score
      @score = geomean assembly.assembly.values.map{ |contig| contig.score }
      return @score
    end

    def assembly_metrics
      @assembly.run unless @assembly.has_run
      @assembly
    end

    def read_metrics left=nil, right=nil
      unless @read_metrics.has_run
        @read_metrics.run(left, right, threads: @threads)
      end
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
