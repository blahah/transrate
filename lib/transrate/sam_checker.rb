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
        diff = position - (length + 1)
        if list[-1][1] == "S"

          list[-1][0] = (list[-1][0].to_i + diff).to_s
          if list[-2][0].to_i > diff
            list[-2][0] = (list[-2][0].to_i - diff).to_s
          elsif list[-2][0].to_i == diff
            list.delete_at(-2) # delete_at changes `list`
          else
          end
        elsif list[-1][1] == "M"
          if list[-1][0].to_i > diff
            list[-1][0] = (list[-1][0].to_i - diff).to_s
            list << [diff.to_s, "S"]
          elsif list[-1][0].to_i == diff
            list[-1][1] = "S"
          end
        end
        cols[5] = list.join("")
        return cols.join("\t")
      else
        return sam
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
              out.write(check(sam1))
              out.write(check(sam))
              @first = true
            end
            @count+=1
          end
        end
      end
    end

  end

end