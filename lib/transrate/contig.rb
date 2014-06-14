require 'forwardable'
require 'inline'

module Transrate

  # A contig an a transcriptome assembly.
  class Contig

    include Enumerable
    extend Forwardable
    def_delegators :@seq, :each, :size, :length

    def init seq
      @seq = seq
    end

    # Length of the contig
    def length
      return @seq.length
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

  end

end
