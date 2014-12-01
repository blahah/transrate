module Transrate

  class SamChecker

    def initialize
      @contigs = {}
      @reference = ""
      @count = 0
      @percent = 0
      @first = true
    end

    def check sam
      cols = sam.split("\t")

      reference = cols[2]
      length = @contigs[reference]

      seq_length = cols[9].length
      position = cols[3].to_i
      cigar = cols[5]
      # this generates a list of pairs in the form [ ["10", "M"], ["1", "D"] ]
      list = cigar.split(/[MDIS]/).zip(cigar.scan(/[MDIS]/))
      list.each_with_index do |a, i|
        c=a[0].to_i
        t=a[1]
        if t=="M" or t=="D"
          position += c
        elsif i==0 and t=="S"
          position += c
        end
      end
      if position > length + 1
        return false
      else
        return true
      end
    end

    def fix_sam input, output
      sam1 = ""
      File.open("#{output}", "wb") do |out|
        File.open("#{input}").each_line do |sam|
          if sam =~ /^@/
            # header
            # @SQ SN:Locus_1_Transcript_13/342_Confidence_1.000_Length_1605 LN:1605
            if sam[0..2]=="@SQ"
              cols = sam.split("\t")
              name = cols[1][3..-1]
              length = cols[2][3..-1].to_i
              @contigs[name] = length
            end
            out.write sam
          else
            # alignment
            if @first
              sam1 = sam.dup
              @first = false
            else
              if check(sam1) and check(sam)
                out.write(sam1)
                out.write(sam)
              end
              @first = true
            end
            @count+=1
          end
        end
      end
    end

  end

end