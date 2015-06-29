require 'set'
require 'crb-blast'

module Transrate

  class ComparativeMetrics

    attr_reader :rbh_per_contig
    attr_reader :rbh_per_reference
    attr_reader :reciprocal_hits
    attr_reader :has_run
    attr_reader :reference_coverage
    attr_reader :comp_stats

    def initialize assembly, reference, threads
      @assembly = assembly
      @reference = reference
      @threads = threads
      @comp_stats = Hash.new
    end

    def run
      crbblast = run_crb_blast
      calculate_reference_coverage crbblast
      @has_run = true
    end

    def calculate_reference_coverage crbblast
      # The reciprocals hash in crb blast has contig names as the key.
      # In order to look up by the reference name we need to reverse this.
      # Scan through the reciprocals and get this Hit objects and add them to
      #   the @reference object for each reference sequence
      get_reference_hits crbblast
      per_query_contig_reference_coverage
      per_target_contig_reference_coverage crbblast
    end

    def get_reference_hits crbblast
      @reference.each do |name, contig|
        contig.hits = []
      end
      crbblast.reciprocals.each do |query_id, list|
        list.each do |hit|
          unless @reference.assembly.key? hit.target
            raise TransrateError.new "#{hit.target} not in reference"
          end
          @reference[hit.target].hits << hit
        end
      end
      @comp_stats[:CRBB_hits] = crbblast.size
      @comp_stats[:n_contigs_with_CRBB] = crbblast.reciprocals.size
      @comp_stats[:p_contigs_with_CRBB] = crbblast.reciprocals.size/@assembly.size.to_f
    end

    def per_query_contig_reference_coverage
      # for each query contig in the @assembly find out how much it covers
      # the reference
      n_refs_with_recip = 0
      total_crbb_hits = 0
      @reference.each do |ref_contig_name, ref_contig|
        ref_contig.hits.each do |hit| # a Hit from query to target
          query_contig_name = hit.query
          unless @assembly.assembly.key? query_contig_name
            raise TransrateError.new "#{query_contig_name} not in assembly"
          end
          @assembly[query_contig_name].has_crb = true
          @assembly[query_contig_name].hits << hit
          raise TransrateError.new "query should not be protein" if hit.qprot
          if hit.tprot
            coverage = 3*hit.alnlen+2 - 3*hit.mismatches - 3*hit.gaps
            coverage /= 3.0*hit.tlen
          else
            coverage = hit.alnlen - hit.mismatches - hit.gaps
            coverage /= hit.tlen.to_f
          end
          @assembly[query_contig_name].reference_coverage = coverage
        end

        if ref_contig.hits.size > 0 # this reference has a crbblast hit
          n_refs_with_recip += 1
        end
        total_crbb_hits += ref_contig.hits.size
      end
      @comp_stats[:rbh_per_reference] = total_crbb_hits / @reference.size.to_f
      @comp_stats[:n_refs_with_CRBB] = n_refs_with_recip
      @comp_stats[:p_refs_with_CRBB] = n_refs_with_recip / @reference.size.to_f
    end

    def per_target_contig_reference_coverage crbblast
      # each target sequence in the reference can have multiple query contigs
      # hit it. to calculate the reference coverage you can't just add up the
      # alignment lengths. you have to make sure that overlaps are taken into
      # account
      coverage_thresholds = [0.25, 0.5, 0.75, 0.85, 0.95]
      coverage_totals = [0, 0, 0, 0, 0]
      prot = crbblast.target_is_prot
      total_coverage = 0
      total_length = 0
      @reference.each do |ref_contig_name, ref_contig|
        if prot
          covered = Array.new(ref_contig.length*3, false)
        else
          covered = Array.new(ref_contig.length, false)
        end
        ref_contig.hits.each_with_index do |hit, i| # a Hit from query to target
          if prot
            if hit.qstart % 3 == 0
              tstart = 3*hit.tstart-4
              tend = 3*hit.tend
            elsif hit.qstart % 3 == 1
              tstart = 3*hit.tstart-2
              tend = 3*hit.tend
            elsif hit.qstart % 3 == 2
              tstart = 3*hit.tstart-3
              tend = 3*hit.tend-1
            end
            if hit.qlen % 3 == 1
              tend += 1
            elsif hit.qlen % 3 == 2
              tend += 2
            end
          else
            tstart = hit.tstart
            tend = hit.tend
          end
          (tstart..tend).each do |b|
            covered[b-1] = true # blast coords are 1 indexed
          end
        end
        coverage = covered.reduce(0) { |sum, v| v ? sum + 1 : sum }
        ref_p = coverage / covered.length.to_f
        ref_contig.reference_coverage = ref_p
        coverage_thresholds.each_with_index do |n, index|
          if ref_p >= n
            coverage_totals[index] += 1
          end
        end

        total_coverage += coverage
        total_length += covered.length
      end

      # calculate proportion of ref sequences with coveragre over thresholds
      coverage_thresholds.each_with_index do |p, i|
        @comp_stats["cov#{(100*p).to_i}".to_sym] = coverage_totals[i]
        @comp_stats["p_cov#{(100*p).to_i}".to_sym] =
                                        coverage_totals[i]/@reference.size.to_f
      end
      @comp_stats[:reference_coverage] = total_coverage / total_length.to_f
    end

    def run_crb_blast
      crbblast = CRB_Blast::CRB_Blast.new @assembly.file, @reference.file
      crbblast.run(1e-5, @threads, true)
      crbblast
    end

  end # ComparativeMetrics

end # Transrate
