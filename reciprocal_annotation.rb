require 'usearch'
require 'rb_hit'

class ReciprocalAnnotation

  attr_reader :l2r_results
  attr_reader :r2l_results
  attr_reader :results

  def initialize assembly, reference
    @assembly = assembly
    @reference = reference
  end

  def run
    self.make_assembly_db
    self.make_reference_db
    left2right, right2left = self.reciprocal_align
    self.parse_results left2right, right2left
    @results
  end

  def make_assembly_db
    unless @assembly.orfs_ublast_db
      assembly_base = File.basename(@assembly.file)
      assembly_orfs = assembly_base + ".orfs"
      @usearch.findorfs @assembly.file, assembly_orfs
      assembly_db = assembly_base + ".udb"
      @usearch.make_udb_ublast assembly_orfs, assembly_db
      @assembly.orfs_ublast_db = assembly_db
    end
  end

  def make_reference_db
    unless @reference.ublast_db
      reference_base = File.basename(@reference.file)
      reference_db = reference_base + ".udb"
      @usearch.make_udb_ublast @reference.file, reference_db
      @reference.ublast_db = reference_db
    end
  end

  def reciprocal_align
    us = Usearch.new
    left2right = us.ublast @assembly.file, @reference.ublast_db
    right2left = us.ublast @reference.file, @assembly.orfs_ublast_db
    [left2right, right2left]
  end

  def parse_results left2right, right2left
    l2r_results = self.load_results_file left2right
    r2l_results = self.load_results_file right2left
    @l2r_hits = self.results_to_hits l2r_results
    @r2l_hits = self.results_to_hits r2l_results
    @results = {}
    @l2r_hits.each_pair do |query, best|
      next if best.nil?
      tbest = @r2l_hits[best.target]
      next if tbest.nil?
      @results[query] = best if query == tbest.target
    end
  end

  def results_to_hits results
    hits = {}
    results.each do |hit|
      if hits.has_key? hit.query
        old_hit = hits[hit.query]
        old_eval, old_bits = old_hit.evalue, old_hit.bitscore
        if hit.bitscore > old_bits
          hits[hit.query] = hit 
        elsif hit.bitscore == old_bits && hit.evalue < old_eval
          hits[hit.query] = hit
        end
      else
        hits[hit.query] = hit
      end
    end
    hits
  end

  def load_results_file file
    results = []
    File.open(file).each_line do |line|
      results << RBHit.new(line.chomp.split("\t"))
    end
    results
  end

end # ReciprocalAnnotation