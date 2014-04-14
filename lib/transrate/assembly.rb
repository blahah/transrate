require 'bio'
require 'bettersam'
require 'csv'
require 'forwardable'

module Transrate

  class Assembly

    include Enumerable
    extend Forwardable
    def_delegators :@assembly, :each, :<<, :size, :length

    attr_accessor :ublast_db
    attr_accessor :orfs_ublast_db
    attr_accessor :protein
    attr_reader :assembly
    attr_reader :has_run

    # number of bases in the assembly
    attr_writer :n_bases

    # assembly filename
    attr_accessor :file

    # assembly n50
    attr_reader :n50

    # Return a new Assembly.
    #
    # - +:file+ - path to the assembly FASTA file
    def initialize file
      @file = file
      @assembly = []
      @n_bases = 0
      Bio::FastaFormat.open(file).each do |entry|
        @n_bases += entry.length
        @assembly << entry
      end
    end

    # Return a new Assembly object by loading sequences
    # from the FASTA-format +:file+
    def self.stats_from20_fasta file
      a = Assembly.new file
      a.basic_stats
    end

    def run threads=8
      stats = self.basic_stats threads
      stats.each_pair do |key, value|
        ivar = "@#{key.gsub(/\ /, '_')}".to_sym
        attr_ivar = "#{key.gsub(/\ /, '_')}".to_sym
        # creates accessors for the variables in stats
        singleton_class.class_eval { attr_accessor attr_ivar }
        self.instance_variable_set(ivar, value)
      end
      @has_run = true
    end

    # Return a hash of statistics about this assembly. Stats are
    # calculated in parallel by splitting the assembly into
    # equal-sized bins and calling Assembly#basic_bin_stat on each
    # bin in a separate thread.
    
    def basic_stats threads=8
      
      # create a work queue to process contigs in parallel
      queue = Queue.new
      
      # split the contigs into equal sized bins, one bin per thread
      binsize = (@assembly.size / threads.to_f).ceil
      Transrate.log.info("Processing #{@assembly.size} contigs in #{threads} bins")
      @assembly.each_slice(binsize) do |bin|
        queue << bin
      end

      # a classic threadpool - an Array of threads that allows
      # us to assign work to each thread and then aggregate their
      # results when they are all finished
      threadpool = []

      # assign one bin of contigs to each thread from the queue.
      # each thread will process its bin of contigs and then wait
      # for the others to finish.
      semaphore = Mutex.new
      stats = []

      threads.times do
        threadpool << Thread.new do |thread|
          # keep looping until we run out of bins
          until queue.empty?

            # use non-blocking pop, so an exception is raised
            # when the queue runs dry
            bin = queue.pop(true) rescue nil
            if bin
              # calculate basic stats for the bin, storing them
              # in the current thread so they can be collected
              # in the main thread.
              bin_stats = basic_bin_stats bin
              semaphore.synchronize { stats << bin_stats }
            end
          end
        end
      end

      # collect the stats calculated in each thread and join
      # the threads to terminate them
      threadpool.each(&:join)

      # merge the collected stats and return then
      merge_basic_stats stats

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
      # x = [90, 70, 50, 30, 10]
      # x2 = x.clone
      # cutoff = x2.pop / 100.0
      # res = []
      n1k = 0
      n10k = 0
      orf_length_sum = 0

      # sort the contigs in ascending length order
      # and iterate over them
      bin.sort_by! { |c| c.seq.size }
      bin.each do |contig|
        
        # increment our long contig counters if this
        # contig is above the thresholds
        n1k += 1 if contig.length > 1_000
        n10k += 1 if contig.length > 10_000

        # add the length of the longest orf to the
        # running total
        orf_length_sum += orf_length(contig.seq)

        # increment the cumulative length and check whether the Nx
        # cutoff has been reached. if it has, store the Nx value and
        # get the next cutoff
        cumulative_length += contig.length
#        if cumulative_length >= @n_bases * cutoff
#          res << contig.length
#          if x2.empty?
#            cutoff=1
#          else
#            cutoff = x2.pop / 100.0
#          end 
#        end
      end

      # calculate and return the statistics as a hash
      mean = cumulative_length / @assembly.size
 #     ns = Hash[x.map { |n| "N#{n}" }.zip(res)]
      {
        "n_seqs" => bin.size,
        "smallest" => bin.first.length,
        "largest" => bin.last.length,
        "n_bases" => n_bases,
        "mean_len" => mean,
        "n_1k" => n1k,
        "n_10k" => n10k,
        "orf_percent" => 300 * orf_length_sum / (@assembly.size * mean)
      }
#      }.merge ns

    end # basic_bin_stats

    def merge_basic_stats stats
      # convert the array of hashes into a hash of arrays
      collect = Hash.new{|h,k| h[k]=[]}
      stats.each_with_object(collect) do |collect, result|
        collect.each{ |k, v| result[k] << v }
      end
      merged = {}
      collect.each_pair do |stat, values|
        if stat == 'orf_percent'  || /N[0-9]{2}/ =~ stat
          # store the mean
          merged[stat] = values.inject(:+) / values.size
        elsif stat == 'smallest'
          merged[stat] = values.min
        elsif stat == 'largest'
          merged[stat] = values.max
        else
          # store the sum
          merged[stat] = values.inject(:+)
        end
      end

      merged
      
    end # merge_basic_stats
     
    # finds longest orf in a sequence
    def orf_length sequence
      longest=0
      (1..6).each do |frame|
        translated = Bio::Sequence::NA.new(sequence).translate(frame)
        translated.split('*').each do |orf|
          if orf.length > longest
            longest=orf.length
          end
        end
      end
      return longest
    end

    # return the number of bases in the assembly, calculating
    # from the assembly if it hasn't already been done.
    def n_bases
      unless @n_bases
        @n_bases = 0
        @assembly.each { |s| @n_bases += s.length }
      end
      @n_bases
    end

    def print_stats
      self.basic_stats.map do |k, v| 
        "#{k}#{" " * (20 - (k.length + v.to_i.to_s.length))}#{v.to_i}"
      end.join("\n")
    end

  end # Assembly

end # Transrate
