require 'bio'

module Transrate

  class ReferenceAlignment

    attr_reader :aligned
    attr_reader :unaligned
    attr_reader :avg_perc_aligned
    attr_reader :genome_stats
    attr_reader :has_run
    attr_reader :blat
    
    def initialize assembly, refgenome, threads, evalue, perc_id_threshold, maxIntron

      @assembly = assembly
      @refgenome = refgenome
      @threads = threads
      @evalue = evalue
      @perc_id_threshold = perc_id_threshold
      @maxIntron = maxIntron
      @genome_stats = Hash.new
    end

    def run
      @blat = blat_results
      @aligned = number_aligned
      @unaligned = number_unaligned
      @avg_perc_aligned = average_percent_aligned
      self.run_genome_stats
      @has_run = true
    end

    def run_genome_stats
      @genome_stats[:aligned] = @aligned
      @genome_stats[:unaligned] = @unaligned
      @genome_stats[:avg_perc_aligned] = @avg_perc_aligned
    end

    def blat_results
      records=[];c=true
      `blat #{@refgenome.file} #{@assembly.file} -out=blast9 -maxIntron=#{@maxIntron} /dev/stdout`.split("\n").each do |line|
        (records << line; c=false) if (c && line[0] != "#")
        c = true if line[0] == "#"
      end
      records=Bio::Blast::Report.new(records.join("\n"))
    end

    def number_aligned
      num = 0
      @blat.each {|x| num += 1 if (x.evalue <= @evalue && x.percent_identity >= @perc_id_threshold)}
      num
    end

    def number_unaligned
      num = 0
      @blat.each {|x| num += 1 if (x.evalue > @evalue || x.percent_identity <= @perc_id_threshold)}
      num
    end

    def average_percent_aligned
      count = 0
      sum = 0
      @blat.each {|x| (count += 1; sum += x.percent_identity) if (x.evalue <= @evalue && x.percent_identity >= @perc_id_threshold)}
      (sum/count).round(2)
    end

  end #ReferenceAlignment
end # Transrate
