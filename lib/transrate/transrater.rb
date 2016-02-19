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
    def initialize(assembly, reference, threads: 1)
      if assembly
        if assembly.is_a?(Assembly)
          @assembly = assembly
        else
          @assembly = Assembly.new(assembly)
        end
        @read_metrics = ReadMetrics.new @assembly
      else
        raise TransrateError.new("assembly is nil")
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

    def classify_contigs cutoff
      @assembly.classify_contigs cutoff
    end

    # Run all analyses
    #
    # @param left [String] path to the left reads
    # @param right [String] path to the right reads
    def run left=nil, right=nil
      assembly_metrics
      if left && right
        read_metrics left, right
        @assembly.classify_contigs @cutoff
      end
      comparative_metrics
    end

    # Reduce all metrics for the assembly to a single quality score
    # by taking the geometric mean of the scores for all contigs
    # and multiplying it by the proportion of fragments whose most likely
    # mapping is consistent with the assembly
    # @return [Integer] the assembly score
    def assembly_score
      if !@score_optimiser
        @score_optimiser = ScoreOptimiser.new(@assembly, @read_metrics)
      end
      return @score_optimiser.raw_score
    end

    def weighted_score
      if !@score_optimiser
        @score_optimiser = ScoreOptimiser.new(@assembly, @read_metrics)
      end
      return @score_optimiser.weighted_score
    end

    def assembly_optimal_score prefix
      if !@score_optimiser
        @score_optimiser = ScoreOptimiser.new(@assembly, @read_metrics)
      end
      return @score_optimiser.optimal_score prefix
    end

    def assembly_metrics
      @assembly.run unless @assembly.has_run
      @assembly
    end

    def read_metrics(left, right)
      unless @read_metrics.has_run
        @read_metrics.run(left, right, threads: @threads)
      end
      if !@score_optimiser
        @score_optimiser = ScoreOptimiser.new(@assembly, @read_metrics)
      end
      @score, @cutoff = @score_optimiser.optimal_score
      @read_metrics
    end

    def good_contigs
      {
        :good_contigs => @assembly.good_contigs,
        :p_good_contigs => @assembly.good_contigs/@assembly.size.to_f
      }
    end

    def comparative_metrics
      @comparative_metrics.run unless @comparative_metrics.has_run
      @comparative_metrics
    end

    def all_metrics left, right
      self.run(left, right)
      all = @assembly.basic_stats
      all.merge!(@read_metrics.read_stats)
      all.merge!(@comparative_metrics.comp_stats)
      all[:score] = @score
      all
    end

  end # Transrater

end # Transrate
