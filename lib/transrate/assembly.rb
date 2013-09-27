require 'bio'
require 'bettersam'
require 'csv'
require 'forwardable'

module Transrate

  class Assembly

    include Enumerable
    extend Forwardable
    def_delegators :@assembly, :each, :<<

    attr_accessor :ublast_db
    attr_accessor :orfs_ublast_db
    attr_accessor :protein

    # number of bases in the assembly
    attr_writer :n_bases

    # assembly filename
    attr_accessor :file

    # assembly n50
    attr_reader :n50

    # Reuturn a new Assembly.
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
      @assembly.sort_by! { |x| x.length }
    end

    # Return a new Assembly object by loading sequences
    # from the FASTA-format +:file+
    def self.stats_from_fasta file
      a = Assembly.new file
      a.basic_stats
    end

    def run
      stats = self.basic_stats
      stats.each_pair do |key, value|
        ivar = "@#{key.gsub(/ /, '_')}".to_sym
        self.instance_variable_set(ivar, value)
      end
    end

    # Return a hash of statistics about this assembly
    def basic_stats
      cumulative_length = 0.0
      # we'll calculate Nx for all these x
      x = [90, 70, 50, 30, 10]
      x2 = x.clone
      cutoff = x2.pop / 100.0
      res = []
      n1k = 0
      n10k = 0
      orf_length_sum = 0
      @assembly.each do |s|
        n1k += 1 if s.length > 1_000
        n10k += 1 if s.length > 10_000
        orf_length_sum += orf_length(s.seq)

        cumulative_length += s.length
        if cumulative_length >= @n_bases * cutoff
          res << s.length
          if x2.empty?
            cutoff=1
          else
            cutoff = x2.pop / 100.0
          end 
        end
      end

      mean = cumulative_length / @assembly.size
      ns = Hash[x.map { |n| "N#{n}" }.zip(res)]
      {
        "n_seqs" => @assembly.size,
        "smallest" => @assembly.first.length,
        "largest" => @assembly.last.length,
        "n_bases" => @n_bases,
        "mean_len" => mean,
        "n_1k" => n1k,
        "n_10k" => n10k,
        "average longest orf" => orf_length_sum/@assembly.size,
        "orf proportion" => (300.0*orf_length_sum)/(@assembly.size*mean.to_f)
      }.merge ns
    end

    # finds longest orf in a sequence
    def orf_length sequence
      longest=0
      (1..6).each do |frame|
        translated = Bio::Sequence::NA.new(sequence).translate(frame)
        translated.split(/\*/).each do |orf|
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
