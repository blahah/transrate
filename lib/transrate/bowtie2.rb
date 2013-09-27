module Transrate

  class Bowtie2

    require 'which'
    include Which

    def initialize
      bowtie2_path = which('bowtie2')
      raise "could not find bowtie2 in the path" if bowtie2_path.empty?
      @bowtie2 = bowtie2_path.first
      bowtie2_build_path = which('bowtie2-build')
      raise "could not find bowtie2-build in the path" if bowtie2_build_path.empty?
      @bowtie2_build = bowtie2_build_path.first
    end

    def map_reads file, left, right=nil, insertsize=200, insertsd=50, outputname=nil
      lbase = File.basename(left)
      rbase = File.basename(right)
      outputname ||= "#{lbase}.#{rbase}.sam"
      realistic_dist = insertsize + (3 * insertsd)
      unless File.exists? outputname
        # construct bowtie command
        bowtiecmd = "#{@bowtie2} -k 3 -p 8 -X #{realistic_dist}" # TODO number of cores should be variable '-p 8'
        bowtiecmd += " --no-unal --local --quiet"
        bowtiecmd += " #{File.basename(file)} -1 #{left}"
        # paired end?
        bowtiecmd += " -2 #{right}" if right
        bowtiecmd += " > #{outputname}"
        # run bowtie
        `#{bowtiecmd}`
      end
      outputname
    end

    def build_index file
      unless File.exists?(file + '.1.bt2')
        `#{@bowtie2_build} --offrate 1 #{file} #{File.basename(file)}`
      end
    end

  end # Bowtie2

end # Transrate