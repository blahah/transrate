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
                  right, insertsize: 200,
                  insertsd: 50, outputname: nil,
                  threads: 8)
      raise Bowtie2Error.new("Index not built") if !@index_built
      lbase = File.basename(left.split(",").first)
      rbase = File.basename(right.split(",").first)
      index = File.basename(@index_name)
      @sam = File.expand_path("#{lbase}.#{rbase}.#{index}.sam")
      realistic_dist = insertsize + (3 * insertsd)
      unless File.exists? @sam
        # construct bowtie command
        bowtiecmd = "#{@bowtie2} --very-sensitive"
        bowtiecmd += " -p #{threads} -X #{realistic_dist}"
        bowtiecmd += " --no-unal"
        bowtiecmd += " --seed 1337"
        bowtiecmd += " -x #{@index_name}"
        bowtiecmd += " -1 #{left}"
        # paired end?
        bowtiecmd += " -2 #{right}" if right
        bowtiecmd += " -S #{@sam}"
        # run bowtie
        runner = Cmd.new bowtiecmd
        runner.run
        # parse bowtie output
        if runner.stderr=~/([0-9]+)\ reads\;\ of\ these\:/
          @read_count = $1.to_i
          # save read_count to file
          File.open("#{@sam}-read_count.txt", "wb") do |out|
            out.write("#{@read_count}\n")
          end
        end
        if !runner.status.success?
          raise Bowtie2Error.new("Bowtie2 failed\n#{runner.stderr}")
        end
      else
        @read_count = 0
        if File.exist?("#{@sam}-read_count.txt")
          @read_count = File.open("#{@sam}-read_count.txt").readlines.join.to_i
        else
          left.split(",").each do |l|
            cmd = "wc -l #{l}"
            count = Cmd.new(cmd)
            count.run
            if count.status.success?
              @read_count += count.stdout.strip.split(/\s+/).first.to_i/4
            else
              logger.warn "couldn't get number of reads from #{l}"
            end
          end
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
