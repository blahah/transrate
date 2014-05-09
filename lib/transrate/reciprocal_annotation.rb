module Transrate

  class ReciprocalAnnotation

    attr_reader :l2r_hits
    attr_reader :r2l_hits
    attr_reader :results

    def initialize assembly, reference
      @assembly = assembly
      @reference = reference
      @usearch = Usearch.new
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
        assembly_dir = File.dirname(@assembly.file)
        assembly_base = File.basename(@assembly.file, ".*")
        assembly_orfs = assembly_base + ".orfs"
        @usearch.findorfs @assembly.file, assembly_orfs
        assembly_db = File.join(assembly_dir, assembly_base + ".udb")
        @usearch.makeudb_ublast assembly_orfs, assembly_db
        @assembly.orfs_ublast_db = assembly_db
      end
    end

    def make_reference_db
      unless @reference.ublast_db
        reference_dir = File.dirname(@reference.file)
        reference_base = File.basename(@reference.file, ".*")
        reference_db = File.join(reference_dir, reference_base + ".udb")
        @usearch.makeudb_ublast @reference.file, reference_db
        @reference.ublast_db = reference_db
        return reference_db
      end
    end

    def reciprocal_align
      left2right = @usearch.ublast @assembly.file, @reference.ublast_db
      right2left = @usearch.ublast @reference.file, @assembly.orfs_ublast_db
      [left2right, right2left]
    end

    def parse_results left2right, right2left
      @l2r_results = self.load_results_file left2right
      @r2l_results = self.load_results_file right2left
      @l2r_hits = self.results_to_hits @l2r_results
      @r2l_hits = self.results_to_hits @r2l_results
      @results = {}
      @l2r_hits.each_pair do |query, best|
        next if best.nil?
        tbest = @r2l_hits[best.target]
        next if tbest.nil?
        if query == tbest.target
        	@results[query] = best 
        end
      end
    end

    # what is this method trying to do? :/
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

end # Transrate
