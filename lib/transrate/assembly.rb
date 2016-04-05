require 'bio'
require 'csv'
require 'forwardable'

module Transrate

  class AssemblyError < TransrateError; end

  # Container for a transcriptome assembly and its associated
  # metadata.
  # @!attribute [r] assembly
  #   @return [Array<Bio::FastaFormat>] the assembly
  # @!attribute [r] has_run
  #   @return [BOOL] whether the basic metrics have been generated
  # @!attribute [w] n_bases
  #   @return [Integer] the number of bases in the assembly
  # @!attribute [rw] file
  #   @return [String] path to the assembly FASTA file
  # @!attribute [r] n50
  #   @return [Integer] assembly n50
  class Assembly

    include Enumerable
    extend Forwardable
    def_delegators :@assembly, :each, :each_value, :<<, :size, :length, :[]

    attr_accessor :file
    attr_reader :assembly
    attr_reader :has_run
    attr_accessor :n_bases
    attr_reader :n50
    attr_accessor :contig_metrics

    # Create a new Assembly.
    #
    # @param file [String] path to the assembly FASTA file
    def initialize file
      @file = File.expand_path file
      unless File.exist? @file
        raise TransrateIOError.new "Assembly file doesn't exist: #{@file}"
      end
      @assembly = {}
      @n_bases = 0
      Bio::FastaFormat.open(file).each do |entry|
        if entry.seq.length == 0
          logger.error "Entry found with no sequence #{entry.entry_id}"
          raise AssemblyError
        end
        @n_bases += entry.length
        contig = Contig.new(entry)
        if @assembly.key?(contig.name)
          logger.error "Non unique fasta identifier found"
          logger.error ">#{contig.name}"
          logger.error "Please make sure there are no duplicate entries in the assembly"
          logger.error "Contig name is taken from before the first | or space"
          logger.error "If you used Trinity, there is a known bug that breaks" +
                       "contig names to make them non-unique."
          logger.error "You can fix your Trinity assembly by replacing | with _"
          logger.error "Example commands to achieve this:"
          logger.error "On OSX:"
          logger.error "`sed 's/\\|/_/' Trinity.fasta > Trinity.fixed.fa`"
          logger.error "On Linux:"
          logger.error "`sed 's_|_-_g' Trinity.fasta > Trinity.fixed.fa`"
          raise AssemblyError
        end
        if contig.name =~ /\,/
          logger.error "Contig names can't contain commas"
          raise AssemblyError
        end
        @assembly[contig.name] = contig
      end
      @contig_metrics = ContigMetrics.new self
    end

    # Generate and store the basic statistics for this assembly
    #
    # @param threads [Integer] number of threads to use
    def run threads=8
      stats = self.basic_stats threads
      stats.each_pair do |key, value|
        ivar = "@#{key.gsub(/\ /, '_')}".to_sym
        attr_ivar = "#{key.gsub(/\ /, '_')}".to_sym
        # creates accessors for the variables in stats
        singleton_class.class_eval { attr_accessor attr_ivar }
        self.instance_variable_set(ivar, value)
      end
      @contig_metrics.run
      @has_run = true
    end

    # Return a hash of statistics about this assembly. Stats are
    # calculated in parallel by splitting the assembly into
    # equal-sized bins and calling Assembly#basic_bin_stat on each
    # bin in a separate thread.
    #
    # @param threads [Integer] number of threads to use
    #
    # @return [Hash] basic statistics about the assembly
    def basic_stats threads=1
      return @basic_stats if @basic_stats
      bin = @assembly.values
      @basic_stats = basic_bin_stats bin
      @basic_stats
    end # basic_stats


    # Calculate basic statistics in an single thread for a bin
    # of contigs.
    #
    # Basic statistics are:
    #
    # - N10, N30, N50, N70, N90
    # - number of contigs >= 1,000 base pairs long
    # - number of contigs >= 10,000 base pairs long
    # - length of the shortest contig
    # - length of the longest contig
    # - number of contigs in the bin
    # - mean contig length
    # - total number of nucleotides in the bin
    # - mean % of contig length covered by the longest ORF
    #
    # @param [Array] bin An array of Bio::Sequence objects
    # representing contigs in the assembly

    def basic_bin_stats bin

      # cumulative length is a float so we can divide it
      # accurately later to get the mean length
      cumulative_length = 0.0

      # we'll calculate Nx for x in [10, 30, 50, 70, 90]
      # to do this we create a stack of the x values and
      # pop the first one to set the first cutoff. when
      # the cutoff is reached we store the nucleotide length and pop
      # the next value to set the next cutoff. we take a copy
      # of the Array so we can use the intact original to collect
      # the results later
      x = [90, 70, 50, 30, 10]
      x2 = x.clone
      cutoff = x2.pop / 100.0
      res = []
      n_under_200, n_over_1k, n_over_10k, n_with_orf, orf_length_sum = 0,0,0,0,0
      # sort the contigs in ascending length order
      # and iterate over them
      bin.sort_by! { |c| c.seq.length }
      bin.each do |contig|
        # increment our long contig counters if this
        # contig is above the thresholds
        if contig.length < 200
          # ignore contigs less than 200 bases,
          # but record how many there are
          n_under_200 += 1
          next
        end
        n_over_1k += 1 if contig.length > 1_000
        n_over_10k += 1 if contig.length > 10_000

        # add the length of the longest orf to the
        # running total
        orf_length = contig.orf_length
        orf_length_sum += orf_length
        # only consider orfs that are realistic length
        # (here we set minimum amino acid length as 50)
        n_with_orf += 1 if orf_length > 149

        # increment the cumulative length and check whether the Nx
        # cutoff has been reached. if it has, store the Nx value and
        # get the next cutoff
        cumulative_length += contig.length
        if cumulative_length >= @n_bases * cutoff
          res << contig.length
          if x2.empty?
            cutoff = 1
          else
            cutoff = x2.pop / 100.0
          end
        end
      end

      # if there aren't enough sequences we might have no value for some
      # of the Nx. Fill the empty ones in with the longest contig length.
      while res.length < x.length do
        res << bin.last.length
      end

      # calculate and return the statistics as a hash
      mean = cumulative_length / @assembly.size
      if @assembly.size * mean == 0
        mean_orf_percent = 0
      else
        mean_orf_percent = 300 * orf_length_sum / (@assembly.size * mean)
      end
      ns = Hash[x.map { |n| "n#{n}" }.zip(res)]
      {
        'n_seqs' => bin.size,
        'smallest' => bin.first.length,
        'largest' => bin.last.length,
        'n_bases' => n_bases,
        'mean_len' => mean,
        'n_under_200' => n_under_200,
        'n_over_1k' => n_over_1k,
        'n_over_10k' => n_over_10k,
        'n_with_orf' => n_with_orf,
        'mean_orf_percent' => mean_orf_percent
      }.merge ns

    end # basic_bin_stats

    def classify_contigs cutoff
      # create hash of file handles for each output
      base = File.basename @file
      files = {}
      %w(good bad).each do |type|
        files[type.to_sym] = File.open("#{type}.#{base}", "wb")
      end
      # loop through contigs writing them out to the appropriate file
      @assembly.each_pair do |name, contig|
        handle = files[contig.classify(cutoff)]
        handle.write contig.to_fasta
      end
      # close all the file handles
      files.each do |type, handle|
        handle.close
      end
      #
      dir = "single_component_bad"
      Dir.mkdir(dir) unless Dir.exist?(dir)
      Dir.chdir(dir) do
        ["cov", "seg","good", "seq"].each do |comp|
          File.open("#{comp}.#{File.basename(@file)}","w") do |out|
            @assembly.each_pair do |name, contig|
              method = "@score_"+comp
              contig.alt_score
              if contig.instance_variable_get(method) > cutoff and
                 contig.score < cutoff
                out.write ">#{name}\n"
                out.write "#{contig.seq.seq}\n"
              end
            end
          end
        end
      end
    end

    def good_contigs
      good = 0
      @assembly.each do |name, contig|
        good += 1 if contig.classification == :good
      end
      good
    end

  end # Assembly

end # Transrate
