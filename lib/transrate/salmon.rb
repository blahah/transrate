module Transrate

  class SalmonError < StandardError
  end

  class Salmon

    def initialize
      which = Cmd.new('which salmon')
      which.run
      if !which.status.success?
        raise SalmonError.new("could not find salmon in the path")
      end
      @salmon = which.stdout.split("\n").first
    end

    def run assembly, bamfile, threads=8
      assembly = assembly.file if assembly.is_a? Assembly
      output = "quant.sf"
      @fin_output = "#{File.basename assembly}_#{output}"

      unless File.exist? @fin_output
        salmon = Cmd.new build_command(assembly, bamfile, threads)
        salmon.run
        unless salmon.status.success?
          logger.warn "salmon failed"
          abort
        end
        File.rename(output, @fin_output)
      end
      return 'postSample.bam'
    end

    def build_command assembly, bamfile, threads=4
      cmd = "#{@salmon} quant"
      cmd << " --libtype IU"
      cmd << " --alignments #{bamfile}"
      cmd << " --targets #{assembly}"
      cmd << " --threads #{threads}"
      cmd << " --useReadCompat"
      cmd << " --useFragLenDist"
      cmd << " --sampleOut"
      cmd << " --sampleUnaligned" # thanks Rob!
      cmd << " --output ."
      cmd
    end

    def load_expression file
      expression = {}
      File.open(file).each do |line|
        if line !~ /^#/
          line = line.chomp.split("\t")
          target = line[0]
          effective_length = line[1]
          effective_count = line[4]
          tpm = line[2]
          expression[target] = {
            :eff_len => effective_length.to_i,
            :eff_count => effective_count.to_f,
            :tpm => tpm.to_f
          }
        end
      end
      expression
    end


  end

end

# # salmon (alignment-based) v0.2.2
# # [ program ] => salmon
# # [ command ] => quant
# # [ libtype ] => { IU }
# # [ alignments ] => { /local_data/cmb211/salmon/SRR037735_1.fastq.SRR037735_2.fastq.transcripts.bam }
# # [ targets ] => { /local_data/cmb211/salmon/transcripts.fa }
# # [ threads ] => { 8 }
# # [ output ] => { rice-salmon }
# # [ sampleOut ] => { }
# # [ useFragLenDist ] => { }
# # [ sampleUnaligned ] => { }
# # Name      Length  TPM       FPKM      NumReads
# scaffold1   1016    549.279   527.364   20690
# scaffold2   1439    598.782   574.892   31945
# scaffold3   783     408.072   391.791   11846
# scaffold4   893     441.382   423.772   14613
# scaffold5   622     494.487   474.758   11403
# scaffold6   2073    4.77214   4.58174   366.764
# scaffold7   1291    4.288     4.11692   205.236
# scaffold8   1355    17.9155   17.2007   900
