
module Transrate

  class ExpressError < StandardError
  end

  class Express

    require 'ostruct'

    attr_reader :fin_output

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
      @fin_output = "#{File.basename assembly}_#{ex_output}"

      unless File.exists? @fin_output
        runner = Cmd.new build_command(assembly, bamfile)
        runner.run
        unless runner.status.success?
          logger.warn "express failed. cleaning sam file and trying again"
          fix_problem_snap_output bamfile
          runner.run
          unless runner.status.success?
            abort "express failed on the cleaned sam file\n#{runner.stderr}"
          end
        end
        File.rename(ex_output, @fin_output)
      end
      return 'hits.1.samp.bam'
    end

    # return the constructed eXpress command
    def build_command assembly, bamfile
      cmd = "#{@express}"
      cmd << " --output-dir ."
      cmd << " --output-align-samp"
      cmd << " --no-update-check"
      cmd << " --additional-online 1"
      cmd << " #{File.expand_path assembly}"
      cmd << " #{File.expand_path bamfile}"
      cmd
    end

    # return a hash of target => effective_count created
    # by parsing the results file
    def load_expression file
      expression = {}
      first = true
      File.open(file).each do |line|
        if first # skip header line
          first = false
          next
        end
        line = line.chomp.split("\t")
        target = line[1]
        effective_length = line[3]
        effective_count = line[7]
        tpm = line[14]
        expression[target] = {
          :eff_len => effective_length.to_i,
          :eff_count => effective_count.to_f,
          :tpm => tpm.to_f
        }
      end
      expression
    end

    def fix_problem_snap_output bam
      # express failed, probably because of temporary snap error
      # convert bam to sam
      sam = "#{File.expand_path(File.basename(bam, File.extname(bam)))}.sam"
      Samtools.run "view -h #{bam} > #{sam}"
      # run sam fixer on sam
      checker = SamChecker.new
      fixed_sam = "#{File.expand_path(File.basename(sam, File.extname(sam)))}.fixed.sam"
      checker.fix_sam(sam, fixed_sam)
      # convert sam to bam
      Samtools.run "view -bS #{fixed_sam} > #{bam}"
      bam
    end

  end # Express

 end # Transrate
