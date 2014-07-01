module Transrate

  class Bowtie2Error < StandardError
  end

  class Bowtie2

    require 'which'
    include Which

    attr_reader :index_name, :sam

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
      lbase = File.basename(left)
      rbase = File.basename(right)
      index = File.basename(@index_name)
      @sam ||= "#{lbase}.#{rbase}.#{index}.sam"
      realistic_dist = insertsize + (3 * insertsd)
      unless File.exists? @sam
        # construct bowtie command
        bowtiecmd = "#{@bowtie2} --very-sensitive"
        bowtiecmd += " -p #{threads} -X #{realistic_dist}"
        bowtiecmd += " --quiet --no-unal"
        bowtiecmd += " --seed 1337"
        bowtiecmd += " -x #{@index_name}"
        bowtiecmd += " -1 #{left}"
        # paired end?
        bowtiecmd += " -2 #{right}" if right
        bowtiecmd += " &> #{@sam}"
        # run bowtie
        `#{bowtiecmd}`
      end
      @sam
    end

    def build_index file
      unless File.exists?(file + '.1.bt2')
        @index_name = File.basename(file).split(".")[0..-2].join(".")
        cmd = "#{@bowtie2_build} --quiet --offrate 1 #{file} #{@index_name}"
        `#{cmd} 2>&1`
        @index_built = true
      end
    end

  end # Bowtie2

end # Transrate
