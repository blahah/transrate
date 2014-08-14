module Transrate

  class Bowtie2Error < StandardError
  end

  class Bowtie2

    require 'which'
    include Which

    attr_reader :index_name, :sam, :read_count

    def initialize
      bowtie2_path = which('bowtie2')
      if bowtie2_path.empty?
        raise Bowtie2Error.new("could not find bowtie2 in the path")
      end
      @bowtie2 = bowtie2_path.first
      bowtie2_build_path = which('bowtie2-build')
      if bowtie2_build_path.empty?
        raise Bowtie2Error.new("could not find bowtie2-build in the path")
      end
      @bowtie2_build = bowtie2_build_path.first
      @index_built = false
      @index_name = ""
    end

    def map_reads(file, left,
                  right, unpaired, library, insertsize: 200,
                  insertsd: 50, outputname: nil,
                  threads: 8)
      raise Bowtie2Error.new("Index not built") if !@index_built
      raise Bowtie2Error.new("Left and right and/or unpaired reads required") unless unpaired || (left && right)
      lbase = File.basename(left.split(",").first) if left
      rbase = File.basename(right.split(",").first) if right
      ubase = File.basename(unpaired.split(",").first) if unpaired
      index = File.basename(@index_name)
      @sam = File.expand_path("#{lbase}.#{rbase}.#{ubase}.#{index}.sam")
      realistic_dist = insertsize + (3 * insertsd)
      unless File.exists? @sam
        # construct bowtie command
        bowtiecmd = "#{@bowtie2} --very-sensitive"
        bowtiecmd += " -p #{threads} -X #{realistic_dist}"
        #bowtiecmd += " --quiet"
        bowtiecmd += " --seed 1337"
        bowtiecmd += " -x #{@index_name}"
        bowtiecmd += " -1 #{left}" if left
        # paired end?
        bowtiecmd += " -2 #{right}" if right
        bowtiecmd += " -U #{unpaired}" if unpaired
		(library == "f" ? bowtiecmd += " --norc" : bowtiecmd += " --norc --#{library}") if library
        bowtiecmd += " -S #{@sam}"
        # run bowtie
        runner = Cmd.new bowtiecmd
        runner.run
        # parse bowtie output
        if runner.stderr=~/([0-9]+)\ reads\;\ of\ these\:/
          @read_count = $1.to_i
        end
        if !runner.status.success?
          raise Bowtie2Error.new("Bowtie2 failed\n#{runner.stderr}")
        end
      else
        cmd = "wc -l #{left}"
        count = Cmd.new(cmd)
        count.run
        if count.status.success?
          @read_count = count.stdout.strip.split(/\s+/).first.to_i/4
        else
          logger.warn "couldn't get number of reads from fastq file"
          @read_count = 0
        end
      end
      @sam
    end

    def build_index file
      @index_name = File.basename(file).split(".")[0..-2].join(".")
      unless File.exists?(@index_name + '.1.bt2')
        cmd = "#{@bowtie2_build} --quiet --offrate 1 #{file} #{@index_name}"
        runner = Cmd.new cmd
        runner.run
        if !runner.status.success?
          msg = "Failed to build Bowtie2 index\n#{runner.stderr}"
          raise Bowtie2Error.new(msg)
        end
      end
      @index_built = true
    end

  end # Bowtie2

end # Transrate
