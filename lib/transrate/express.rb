module Transrate

  require 'which'

  class Express

    # return an Express object
    def initialize
      express_path = Which::which('express')
      raise "could not find eXpress in the path" if express_path.empty?
      @express = express_path.first
    end

    # return hash of expression for each sequenceID
    # in the assembly fastafile
    def quantify_expression assembly, samfile
      assembly = assembly.file if assembly.is_a? Assembly
      cmd = "#{@express} --no-bias-correct #{File.expand_path assembly} "
      cmd << " #{File.expand_path samfile}"
      ex_output = 'results.xprs'
      fin_output = "#{assembly}_#{ex_output}"
      unless File.exists? fin_output
        stdout = `#{cmd} 2>&1`.split(/\n/)[1..30].join("\n")
        File.rename(ex_output, fin_output)
      end
      expression = {}
      File.open(fin_output).each do |line|
        line = line.chomp.split("\t")
        target = line[1]
        effective_count = line[7]
        expression[target] = effective_count
      end
      expression
    end

  end # Express

 end # Transrate
