
module Transrate

  class ExpressError < StandardError
  end

  class Express

    # return an Express object
    def initialize
      which = Cmd.new('which express')
      which.run
      if !which.status.success?
        raise ExpressError.new("could not find express in the path")
      end
      @express = which.stdout.split("\n").first
    end

    # return struct containing:
    #   results_file => path to the express results TSV
    #   expression   => a hash of target => effective_count
    #   align_samp   => path to the sampled alignments file
    def run assembly, bamfile
      assembly = assembly.file if assembly.is_a? Assembly

      ex_output = 'results.xprs'
      fin_output = "#{File.basename assembly}_#{ex_output}"

      unless File.exists? fin_output
        runner = Cmd.new build_command(assembly, bamfile)
        runner.run
        unless runner.status.success?
          raise ExpressError.new("Express failed\n#{runner.stderr}")
        end
        File.rename(ex_output, fin_output)
      end

      OpenStruct.new(:results_file => fin_output,
                     :expression => load_expression(fin_output),
                     :align_samp => 'hits.1.samp.bam')
    end

    # return the constructed eXpress command
    def build_command assembly, bamfile
      cmd = "#{@express}"
      cmd << " #{File.expand_path assembly}"
      cmd << " #{File.expand_path bamfile}"
      cmd << " --output-dir ."
      cmd << " --output-align-samp"
      cmd << " --no-update-check"
      cmd << " --additional-online 1"
      cmd
    end

    # return a hash of target => effective_count created
    # by parsing the results file
    def load_expression file
      expression = {}
      first = true
      File.open(file).each do |line|
        if first
          first = false
          next
        end
        line = line.chomp.split("\t")
        target = line[1]
        effective_count = line[7]
        expression[target] = effective_count.to_f
      end
      expression
    end

  end # Express

 end # Transrate
