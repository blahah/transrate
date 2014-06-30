module Transrate

  class Bowtie2Error < StandardError
  end

  class Bowtie2

    require 'which'
    include Which

    attr_reader :index_name, :sam

    def initialize
      bowtie2_path = which('bowtie2')
      raise Bowtie2Error.new("could not find bowtie2 in the path") if bowtie2_path.empty?
      @bowtie2 = bowtie2_path.first
      bowtie2_build_path = which('bowtie2-build')
      raise Bowtie2Error.new("could not find bowtie2-build in the path") if bowtie2_build_path.empty?
      @bowtie2_build = bowtie2_build_path.first
      @index_built = false
      @index_name = ""
    end

    def map_reads file, left, right=nil, insertsize=200, insertsd=50, sam=nil
      raise Bowtie2Error.new("Index not built") if !@index_built
      lbase = File.basename(left)
      rbase = File.basename(right)
      index = File.basename(@index_name)
      dir = File.dirname(left)
      @sam ||= "#{dir}/#{lbase}.#{rbase}.#{index}.sam"
      realistic_dist = insertsize + (3 * insertsd)
      unless File.exists? @sam
        # construct bowtie command
        bowtiecmd = "#{@bowtie2} --very-sensitive-local -k 10 -p 8 "
        bowtiecmd += " -X #{realistic_dist}"
        bowtiecmd += " --no-unal --quiet"
        bowtiecmd += " -x #{@index_name} -1 #{left}"
        # paired end?
        bowtiecmd += " -2 #{right}" if right
        bowtiecmd += " > #{@sam}"
        # run bowtie
        `#{bowtiecmd}`
      end
      @sam
    end

    def build_index file
      unless File.exists?(file + '.1.bt2')
        @index_name = file.split(".")[0..-2].join(".")
        cmd = "#{@bowtie2_build} --quiet --offrate 1 #{file} #{@index_name}"
        `#{cmd}`
        @index_built = true
      end
    end

  end # Bowtie2

end # Transrate
