require 'forwardable'

module Transrate

  # A contig in a transcriptome assembly.
  class Contig

    include Enumerable
    extend Forwardable
    def_delegators :@seq, :size, :length
    attr_accessor :seq, :name
    # read-based metrics
    attr_accessor :eff_length, :eff_count, :tpm
    attr_accessor :coverage, :uncovered_bases, :p_uncovered_bases
    attr_accessor :p_seq_true
    attr_accessor :low_uniqueness_bases, :in_bridges
    attr_accessor :p_good, :p_not_segmented, :good, :classification
    # reference-based metrics
    attr_accessor :has_crb, :reference_coverage
    attr_accessor :hits
    attr_accessor :score_cov, :score_seg, :score_good, :score_seq

    def initialize(seq, name: nil)
      # fix null bytes in the nucleotide sequence
      seq.seq.gsub!("\0", "")
      # trim trailing semicolons (because BLAST strips them)
      if seq.respond_to?(:entry_id)
        seq.entry_id.gsub!(/;$/, '')
      end
      @seq = seq
      @seq.data = nil # no need to store raw fasta string
      @name = seq.respond_to?(:entry_id) ? seq.entry_id : name
      @hits = []
      @reference_coverage = 0
      @has_crb = false
      @in_bridges = 0
      @p_good = 0
      @p_seq_true = 0
      @uncovered_bases = length
      @p_uncovered_bases = 1
      @p_not_segmented = 1
      @score = -1
      @good = 0
      @coverage = 0
      @classification = :unknown
    end

    def each &block
      @seq.seq.each_char &block
    end

    # Get all metrics available for this contig
    def basic_metrics
      basic = {
        :length => length,
        :prop_gc => prop_gc,
        :orf_length => orf_length
      }
    end

    def read_metrics
      {
        :in_bridges => in_bridges,
        :p_good => p_good,
        :p_bases_covered => p_bases_covered,
        :p_seq_true => p_seq_true,
        :score => score,
        :p_not_segmented => p_not_segmented,
        :eff_length => eff_length,
        :eff_count => eff_count,
        :tpm => tpm,
        :coverage => coverage,
        :sCnuc => p_seq_true,
        :sCcov => p_bases_covered,
        :sCord => p_good,
        :sCseg => p_not_segmented
      }
    end

    def comparative_metrics
      reference = @has_crb ? {
        :has_crb => has_crb,
        :reference_coverage => reference_coverage,
        :hits => hits.map{ |h| h.target }.join(";")
      } : {
        :has_crb => false,
        :reference_coverage => "NA",
        :hits => "NA"
      }
    end

    # Base composition of the contig
    #
    # If called and the instance variable @base_composition is nil
    # then call the c method to count the bases and dibases in the sequence
    # then get the info out of the c array and store it in the hash
    # then if it is called again just return the hash as before
    def base_composition
      if @base_composition
        return @base_composition
      end
      # else run the C method
      composition(@seq.seq)
      alphabet = ['a', 'c', 'g', 't', 'n']
      @base_composition = {}
      @dibase_composition = {}
      bases = []
      dibases = []
      alphabet.each do |c|
        bases << "#{c}".to_sym
      end
      alphabet.each do |c|
        alphabet.each do |d|
          dibases << "#{c}#{d}".to_sym
        end
      end
      bases.each_with_index do |a,i|
        @base_composition[a] = base_count(i)
      end
      dibases.each_with_index do |a,i|
        @dibase_composition[a] = dibase_count(i)
      end
      return @base_composition
    end

    # Dibase composition of the contig
    def dibase_composition
      if @dibase_composition
        return @dibase_composition
      end
      base_composition
      @dibase_composition
    end

    # Number of bases that are C
    def bases_c
      base_composition[:c]
    end

    # Proportion of bases that are C
    def prop_c
      bases_c / length.to_f
    end

    # Number of bases that are G
    def bases_g
      base_composition[:g]
    end

    # Proportion of bases that are G
    def prop_g
      bases_g / length.to_f
    end

    # Number of bases that are A
    def bases_a
      base_composition[:a]
    end

    # Proportion of bases that are A
    def prop_a
      bases_a / length.to_f
    end

    # Number of bases that are T
    def bases_t
      base_composition[:t]
    end

    # Proportion of bases that are T
    def prop_t
      bases_t / length.to_f
    end

    def bases_n
      base_composition[:n]
    end

    def prop_n
      bases_n / length.to_f
    end

    # GC
    def bases_gc
      bases_g + bases_c
    end

    def prop_gc
      prop_g + prop_c
    end

    # Find the longest orf in the contig
    def orf_length
      return @orf_length if @orf_length
      @orf_length = longest_orf(@seq.seq) # call to C
      return @orf_length
    end

    def p_bases_covered
      1 - p_uncovered_bases
    end

    def uncovered_bases= n
      @uncovered_bases = n
      @p_uncovered_bases = n / length.to_f
    end

    # Contig score (product of all score components)
    def score
      return @score if @score != -1
      prod =
        [p_bases_covered, 0.01].max.to_f * # proportion of bases covered
        [p_not_segmented, 0.01].max.to_f * # prob contig has 0 changepoints
        [p_good, 0.01].max.to_f * # proportion of reads that mapped good
        [p_seq_true, 0.01].max.to_f # scaled 1 - mean per-base edit distance
      @score = [prod, 0.01].max
    end

    def alt_score
      hash = {
        :cov =>  [p_bases_covered, 0.01].max.to_f,
        :seg =>  [p_not_segmented, 0.01].max.to_f,
        :good => [p_good, 0.01].max.to_f,
        :seq =>  [p_seq_true, 0.01].max.to_f
      }
      hash.keys.each do |miss|
        prod = 1
        hash.each do |key,score|
          if key!=miss
            prod *= score
          end
        end
        score_part = ("score_" + miss.to_s).to_sym
        instance_variable_set("@#{score_part}", prod)
      end
    end

    # Classify the contig into one of the following classes:
    def classify cutoff
      if score >= cutoff
        @classification = :good
      else
        @classification = :bad
      end
      return @classification
    end

    def to_fasta
      @seq.seq.to_fasta(@name)
    end

  end

end
