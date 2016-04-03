module Transrate

  class SalmonError < TransrateError
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
      sampled_bam = "postSample.bam"
      @fin_output = "#{File.basename assembly}_#{output}"
      unless File.exist? @fin_output
        salmon = Cmd.new build_command(assembly, bamfile, threads)
        salmon.run
        unless salmon.status.success?
          logger.error salmon.stderr
          raise SalmonError.new("Salmon failed")
        end
        unless File.exist?(sampled_bam)
          logger.error salmon.stderr
          raise SalmonError.new("#{sampled_bam} not created")
        end
        File.rename(output, @fin_output)
      end
      return sampled_bam
    end

    def build_command assembly, bamfile, threads=4
      cmd = "#{@salmon} quant"
      cmd << " --libType IU"
      cmd << " --alignments #{bamfile}"
      cmd << " --targets #{assembly}"
      cmd << " --threads #{threads}"
      cmd << " --sampleOut"
      cmd << " --sampleUnaligned" # thanks Rob!
      cmd << " --output ."
      cmd << " --useErrorModel"
      cmd << " --biasCorrect"
      cmd << " --noEffectiveLengthCorrection"
      cmd << " --useFSPD"
      cmd
    end

    def load_expression file
      expression = {}
      first = true
      File.open(file).each do |line|
        if first
          first = false
          next
        end
        line = line.chomp.split("\t")
        unless line.length == 5
          raise SalmonError.new("Salmon output file should have 5 columns " +
            "but it had #{line.length}\n" +
            "Please check you are using the correct version of Salmon")
        end
        target = line[0]
        # we ignore length in the second column - we already know it
        effective_length = line[2]
        tpm = line[3]
        effective_count = line[4]

        expression[target] = {
          :eff_len => effective_length.to_i,
          :eff_count => effective_count.to_f,
          :tpm => tpm.to_f
        }
      end
      expression
    end


  end

end
