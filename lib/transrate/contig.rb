require 'forwardable'
require 'inline'

module Transrate

  # A contig in a transcriptome assembly.
  class Contig

    include Enumerable
    extend Forwardable
    def_delegators :@seq, :size, :length
    attr_accessor :seq, :name
    # read-based metrics
    attr_accessor :coverage, :uncovered_bases, :mean_coverage
    # reference-based metrics
    attr_accessor :has_crb, :is_chimera, :collapse_factor, :reference_coverage
    attr_accessor :hits

    def initialize(seq, name: nil)
      @seq = seq
      @name = seq.respond_to?(:entry_id) ? seq.entry_id : name
      @hits = []
      @reference_coverage = 0
      @collapse_factor = 0
      @is_chimera = false
    end

    def each &block
      @seq.seq.each_char &block
    end

    # Get all metrics available for this contig
    def metrics
      basic = {
        :length => length,
        :prop_gc => prop_gc,
        :gc_skew => gc_skew,
        :at_skew => at_skew,
        :cpg_count => cpg_count,
        :cpg_ratio => cpg_ratio,
        :orf_length => orf_length,
        :linguistic_complexity_6 => linguistic_complexity(6)
      }
      read = @coverage ? {
        :uncovered_bases => uncovered_bases,
        :mean_coverage => mean_coverage,
        :in_bridge => in_bridge
      } : {}
      reference = @comp ? {
        :has_crb => has_crb,
        :collapse_factor => collapse_factor,
        :reference_coverage => reference_coverage,
        :is_chimera => is_chimera,
        :hits => hits.map{ |h|  }
      } : {}
      basic.merge(read).merge(reference)
    end

    # Base composition of the contig
    def base_composition
      if @base_composition
        return @base_composition
      end
      base_comp = {
        :a => 0,
        :t => 0,
        :c => 0,
        :g => 0,
        :n => 0
      }
      dibase_comp = {
        :cg => 0
      }
      last_base = nil
      @seq.seq.each_char do |base|
        # single bases
        key = base.downcase.to_sym
        base_comp[key] += 1
        if last_base
          # pairs of bases
          dikey = "#{last_base}#{base}".downcase.to_sym
          if dibase_comp[dikey]
            dibase_comp[dikey] += 1
          else
            dibase_comp[dikey] = 1
          end
        end
        last_base = base
      end
      @base_composition = base_comp
      @dibase_composition = dibase_comp
      return base_comp
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
      (bases_g - bases_c) / (bases_g + bases_c).to_f
    end

    # AT skew
    def at_skew
      (bases_a - bases_t) / (bases_a + bases_t).to_f
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
      longest = longest_orf @seq.seq
      return longest
    end

    # Inlined C longest-ORF function
    inline do |builder|
      builder.c <<SRC
        static
        void
        longest_orf(VALUE _s) {
          int i,sl,longest=0;
          int len[6];
          char * c_str;

          sl = RSTRING_LEN(_s);
          c_str = StringValueCStr(_s);
          for (i=0;i<6;i++) {
            len[i]=0;
          }
          for (i=0;i<sl-2;i++) {
            if (c_str[i]=='T' &&
              ((c_str[i+1]=='A' && c_str[i+2]=='G') ||
              (c_str[i+1]=='A' && c_str[i+2]=='A') ||
              (c_str[i+1]=='G' && c_str[i+2]=='A'))) {
              if (len[i%3] > longest) {
                longest = len[i%3];
              }
              len[i%3]=0;
            } else {
              len[i%3]++;
            }
            if (c_str[i+2]=='A' &&
              ((c_str[i]=='C' && c_str[i+1]=='T') ||
              (c_str[i]=='T' && c_str[i+1]=='T') ||
              (c_str[i]=='T' && c_str[i+1]=='C'))) {
              if (len[3+i%3] > longest) {
                longest = len[3+i%3];
              }
              len[3+i%3]=0;
            } else {
              len[3+i%3]++;
            }
          }
          if (len[i%3] > longest) {
            longest = len[i%3];
          }
          if (len[3+i%3] > longest) {
            longest = len[3+i%3];
          }
          return INT2NUM(longest);
        }
SRC
    end

    inline do |builder|
      builder.c <<SRC
        VALUE
        kmer_count(VALUE _k, VALUE _s) {
          int n, i, start, k, len, h, size = 0;
          char * c_str;
          len = RSTRING_LEN(_s);
          c_str = StringValueCStr(_s);
          k = NUM2INT(_k);
          size = 1;
          for(h=0;h<k;h++) {
            size *= 4;
          }
          short set[size];
          for(start=0;start<size;start++) {
            set[start]=0;
          }
          for(start=0; start<len-k+1; start++) {
            i = 0;
            h = 0;
            n = 0;
            for(i = start; i < start+k; i++) {
              if (c_str[i]=='N') {
                n++;
              }
              if (c_str[i]=='A') {
                h = h << 2;
                h += 0;
              }
              if (c_str[i]=='C') {
                h = h << 2;
                h += 1;
              }
              if (c_str[i]=='G') {
                h = h << 2;
                h += 2;
              }
              if (c_str[i]=='T') {
                h = h << 2;
                h += 3;
              }
            }
            if (n==0) {
              set[h] += 1;
            }
          }
          i = 0; // count how many in array are set //
          for(start = 0; start < size; start++) {
            if (set[start]>0) {
              i++;
            }
          }
          return INT2NUM(i);
        }
SRC
    end

    def linguistic_complexity k
      return kmer_count(k, @seq.seq)/(4**k).to_f
    end
  end

end
