require 'forwardable'
require 'inline'

module Transrate

  # A contig in a transcriptome assembly.
  class Contig

    include Enumerable
    extend Forwardable
    def_delegators :@seq, :size, :length
    attr_accessor :seq, :name, :coverage

    def initialize(seq, name: nil)
      @seq = seq
      @name = seq.respond_to?(:entry_id) ? seq.entry_id : name
    end

    def each &block
      @seq.seq.each_char &block
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
      @dibase_composition={}
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

    # GC skew
    def gc_skew
      prop_gc / (prop_a + prop_t + prop_gc)
    end

    # AT skew
    def at_skew
      prop_a + prop_t / (prop_a + prop_t + prop_gc)
    end

    # CpG count
    def cpg_count
      dibase_composition[:cg]
    end

    # CpG (C-phosphate-G) ratio
    def cpg_ratio
      dibase_composition[:cg] / (prop_c * prop_g)
    end

    # Find the longest orf in the contig
    def orf_length
      return longest_orf(@seq.seq) # call to C
    end

    def linguistic_complexity k
      return kmer_count(k, @seq.seq)/(4**k).to_f # call to C
    end
  end

end
