module Transrate

  class Bowtie2Error < StandardError
  end

  class Bowtie2

    attr_reader :index_name, :sam, :read_count

    def initialize
      @bowtie2 = get_bin_path("bowtie2")
      @bowtie2_build = get_bin_path("bowtie2-build")

      @index_built = false
      @index_name = ""
    end

    def get_bin_path bin
      which_bin = Cmd.new("which #{bin}")
      which_bin.run
      if !which_bin.status.success?
        raise IOError.new("ReadMetrics: could not find #{bin} in path")
      end
      which_bin.stdout.split("\n").first
    end

    def map_reads(file, left, right,
                  insertsize: 200, insertsd: 50, threads: 8)
      raise Bowtie2Error.new("Index not built") if !@index_built
      lbase = File.basename(left.split(",").first)
      rbase = File.basename(right.split(",").first)
      index = File.basename(@index_name)
      @sam = File.expand_path("#{lbase}.#{rbase}.#{index}.sam")
      realistic_dist = insertsize + (5 * insertsd)
      unless File.exists? @sam
        # construct bowtie command
        bowtiecmd = "#{@bowtie2} --very-sensitive"
        bowtiecmd << " -x #{@index_name}"
        bowtiecmd << " -a"
        bowtiecmd << " -p #{threads} -X #{realistic_dist}"
        bowtiecmd << " --no-unal"
        bowtiecmd << " --seed 1337"
        bowtiecmd << " --rdg 6,5"
        bowtiecmd << " --rfg 6,5"
        bowtiecmd << " --score-min L,-.6,-.4"
        bowtiecmd << " -1 #{left}"
        # paired end?
        bowtiecmd << " -2 #{right}" if right
        bowtiecmd << " -S #{@sam}"
        bowtiecmd << " --reorder"
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
              File.open("#{@sam}-read_count.txt", "wb") do |out|
                out.write("#{@read_count}\n")
              end
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