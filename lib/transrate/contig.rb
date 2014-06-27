require 'forwardable'
require 'inline'

module Transrate

  # A contig an a transcriptome assembly.
  class Contig

    include Enumerable
    extend Forwardable
    def_delegators :@seq, :each, :size, :length
    attr_accessor :bases_c, :bases_g, :bases_a, :bases_t,
                  :prop_c, :prop_g, :prop_a, :prop_t,
                  :gc, :gc_skew, :at_skew, :cpg,
                  :linguistic_complexity, :zip_size,
                  :bases_n, :prop_n

    def init seq
      @seq = seq
    end

    # Base composition of the contig
    def base_composition
      if @base_composition
        return @base_composition
      end
      base_comp = {}
      dibase_comp = {}
      last_base = nil
      self.each do |base|
        key = base.downcase.to_sym
        base_comp[dikey] ||= 1
        base_comp[dikey] += 1
        if last_base
          dikey = "#{last_base}#{base}".downcase.to_sym
          dibase_comp[dikey] ||= 1
          dibase_comp[dikey] += 1
        end
        last_base = base
      end
      @base_composition = base_comp
      @dibase_composition = dibase_comp
    end

    # Dibase composition of the contif
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

    # GC
    def gc
      prop_g + prop_c
    end

    # GC skew
    def gc_skew
      gc / (prop_a + prop_t + gc)
    end

    # AT skew
    def at_skew
      prop_a + prop_t / (prop_a + prop_t + gc)
    end

    # CpG (C-phosphate-G) ratio
    def cpg_ratio
      @dibase_composition[:cg] / (bases_c + bases_g) * length
    end

    # Find the longest orf in the contig
    def orf_length sequence
      longest = longest_orf(sequence)
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

    def gc_skew

    end

    def at_skew

    end

    def cpg

    end

    def bases_n

    end

    def proportion_n

    end

    def linguistic_complexity

    end

  end

end
