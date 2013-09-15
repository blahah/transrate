#! /usr/bin/env ruby

require 'bio'
require 'better_sam'
require 'csv'
require 'trollop'

class Assembly

  # number of bases in the assembly
  attr_writer :n_bases

  # assembly filename

  # Reuturn a new Assembly.
  #
  # - +:assembly+ - an array of Bio::Sequences
  def initialize file
    @file = File.expand_path file
    @assembly = []
    @n_bases = 0
    Bio::FastaFormat.open(file).each do |entry|
      @n_bases += entry.length
      @assembly << entry.to_seq
    end
    @assembly.sort_by! { |x| x.length }
  end

  # Return a new Assembly object by loading sequences
  # from the FASTA-format +:file+
  def self.stats_from_fasta file
    a = Assembly.new file
    a.basic_stats
  end

  # Return a hash of statistics about this assembly
  def basic_stats
    cum = 0.0
    x = [90, 70, 50, 30, 10]
    x2 = x.clone
    cutoff = x2.pop / 100.0
    res = []
    n1k = 0
    n10k = 0
    @assembly.each do |s|
      newcum = cum + s.length
      prop = newcum / self.n_bases
      n1k += 1 if s.length > 1_000
      n10k += 1 if s.length > 10_000
      if prop >= cutoff
        res << s.length
        break if x2.empty?
        cutoff = x2.pop / 100.0
      end
      cum = newcum
    end
    mean = cum / @assembly.size
    ns = Hash[x.map { |n| "N#{n}" }.zip(res)]
    {
      "n_seqs" => @assembly.size,
      "smallest" => @assembly.first.length,
      "largest" => @assembly.last.length,
      "n_bases" => @n_bases,
      "mean_len" => mean,
      "n > 1k" => n1k,
      "n > 10k" => n10k
    }.merge ns
  end

  # return the number of bases in the assembly, calculating
  # from the assembly if it hasn't already been done.
  def n_bases
    unless @n_bases
      @n_bases = 0
      @assembly.each { |s| @n_bases += s.length }
    end
    @n_bases
  end

  def map_reads left, right=nil, insertsize=200, insertsd=50, outputname=nil
    self.build_bowtie_index
    lbase = File.basename(left)
    rbase = File.basename(right)
    outputname ||= "#{lbase}.#{rbase}.sam"
    realistic_dist = insertsize + (3 * insertsd)
    unless File.exists? joinedname
      # construct bowtie command
      bowtiecmd = "bowtie2 -k 3 -p 8 -X #{realistic_dist} --no-unal --local --quiet #{@filename} -1 #{left}"
      # paired end?
      bowtiecmd += " -2 #{right}" if right
      # other functions may want the output, so we save it to file
      bowtiecmd += " > #{outputname}"
      # run bowtie
      `#{bowtiecmd}`
    end
  end

  def build_bowtie_index
    unless File.exists?(@assembly + '.1.bt2')
      `bowtie2-build --offrate 1 #{@assembly} #{@assembly_name}`
    end
  end

  def analyse_read_mappings filename, insertsize=200, insertsd=50, bridge=true
    diagnostics = {
      :total => 0,
      :good => 0,
      :bad => 0,
      :both_mapped => 0,
      :not_both_mapped => 0,
      :proper_pair => 0,
      :improper_pair => 0,
      :proper_orientation => 0,
      :improper_orientation => 0,
      :same_contig => 0,
      :realistic => 0,
      :unrealistic => 0
    }
    bridges = {}
    realistic_dist = insertsize + (3 * insertsd)
    if File.exists?(filename) && File.size(filename) > 0
      ls = BetterSam.new
      rs = BetterSam.new
      sam = File.open(filename).each_line
      sam.each_slice(2) do |l, r|
        next if l[0] == "@"
        if l && ls.parse_line(l) # Returns false if line starts with @ (a header line)
          if r && rs.parse_line(r)
            next unless ls.primary_aln?
            diagnostics[:total] += 1
            if ls.read_unmapped? || rs.read_unmapped?
              diagnostics[:not_both_mapped] += 1
            else
              # reads are paired
              diagnostics[:both_mapped] += 1
              if ls.read_properly_paired?
                # mapped in proper pair
                diagnostics[:proper_pair] += 1
                if (!ls.read_reverse_strand? && ls.mate_reverse_strand? || 
                 (ls.read_reverse_strand? && !ls.mate_reverse_strand?))
                  # mates in proper orientation
                  diagnostics[:proper_orientation] += 1
                  diagnostics[:good] += 1
                else
                  # mates in wrong orientation
                  diagnostics[:improper_orientation] += 1
                  diagnostics[:bad] += 1
                end
              else
                # not mapped in proper pair
                diagnostics[:improper_pair] += 1
                # both read and mate are mapped
                diagnostics[:both_mapped] += 1
                if ls.chrom == rs.chrom
                  # both on same contig
                  diagnostics[:same_contig] += 1
                  begin
                    if Math.sqrt((ls.pos - rs.pos) ** 2) < ls.seq.length
                      # overlap is realistic
                      diagnostics[:realistic] += 1
                      if (!ls.read_reverse_strand? && ls.mate_reverse_strand? || 
                       (ls.read_reverse_strand? && !ls.mate_reverse_strand?))
                        # mates in proper orientation
                        diagnostics[:proper_orientation]
                        diagnostics[:good] += 1
                      else
                        # mates in wrong orientation
                        diagnostics[:improper_orientation]
                        diagnostics[:bad] += 1
                      end
                    else
                      # overlap not realistic
                      diagnostics[:unrealistic] += 1
                      diagnostics[:bad] += 1
                    end
                  rescue
                    puts ls.pos
                    puts rs.pos
                  end
                else
                  # mates on different contigs
                  # are the mapping positions within a realistic distance of
                  # the ends of contigs?
                  ldist = [ls.pos, ls.seq.length - ls.pos].min
                  rdist = [rs.pos, rs.seq.length - rs.pos].min
                  if ldist + rdist <= realistic_dist
                    # increase the evidence for this bridge
                    key = [ls.chrom, rs.chrom].sort.join("<>").to_sym
                    if bridges.has_key? key
                      bridges[key] += 1
                    else
                      bridges[key] = 1
                    end
                    diagnostics[:realistic] += 1
                    diagnostics[:good] += 1
                  else
                    diagnostics[:unrealistic] += 1
                    diagnostics[:bad] += 1
                  end
                end
              end
            end
          end
        end
      end
      diagnostics[:supported_bridges] = 0
      if bridge
        CSV.open('supported_bridges.csv', 'w') do |f|
          bridges.each_pair do |b, count|
            start, finish = b.to_s.split('<>')
            if count > 1
              f << [start, finish, count]
              diagnostics[:supported_bridges] += 1
            end
          end
        end
      end
      puts "#{bridges.size} bridges with 2 or more supporting read pairs written to supported_bridges.csv"
      diagnostics[:result] = diagnostics[:bad]
    end
    return diagnostics
  end

  def print_stats
    self.basic_stats.map do |k, v| 
      "#{k}#{" " * (20 - (k.length + v.to_i.to_s.length))}#{v.to_i}"
    end.join("\n")
  end

end # Assembly

def pretty_print_hash hash, width
  hash.map{ |k, v| "#{k.to_s}#{" " * (width - (k.length + v.to_i.to_s.length))}#{v.to_i}" }.join("\n")
end

puts "|" *  40
puts "\n"
a = Assembly.new ARGV.first
puts "Basic assembly stats:"
puts "-" *  40
s = a.basic_stats
puts pretty_print_hash s, 40
puts "|" *  40
puts "\n"
puts "Read mapping statistics:"
puts "-" *  40
brm = a.analyse_read_mappings ARGV[1]
puts pretty_print_hash brm, 40

