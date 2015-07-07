require 'helper'
require 'crb-blast'

module Transrate
  class ComparativeMetrics
    attr_reader :assembly
    attr_reader :reference
    attr_reader :crbblast
  end
end

module CRB_Blast
  class CRB_Blast
    def add_missing
      @missed.each do |query_id, missed|
        missed.each do |hit|
          @reciprocals[hit.query] ||= []
          @reciprocals[hit.query] << hit
        end
      end
    end
  end
end

class Tester
  def self.testpath file
    return File.join(File.dirname(__FILE__), 'data', file)
  end

  def self.run_comp_metrics(query, target)
    querypath = testpath(query)
    targetpath = testpath(target)
    @assembly = Transrate::Assembly.new(querypath)
    @reference = Transrate::Assembly.new(targetpath)
    @comp = Transrate::ComparativeMetrics.new(@assembly, @reference, 1)
    @comp.run
    return @comp
  end
end

class TestCompMetrics2 < MiniTest::Test


  context "ComparativeMetrics" do

    setup do
    end

    should "01 should run" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          comp = Tester.run_comp_metrics("test_contig_nc1.fa", "test_reference_nc1.fa")
          assert comp.has_run
        end
      end
    end

    should "01-1n should get reference hits" do
      # The reciprocals hash in crb blast has contig names as the key.
      # In order to look up by the reference name we need to reverse this.
      # Scan through the reciprocals and get this Hit objects and add them to
      #   the @reference object for each reference sequence
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("test_contig_nc1.fa")
          target = Tester.testpath("test_reference_nc1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          assert_equal 1, comp.reference["reference1"].hits.size, "size of reference hits list"
          assert_equal "contig1", comp.reference["reference1"].hits[0].query
          assert_equal "reference1", comp.reference["reference1"].hits[0].target
        end
      end
    end

    should "01-1n get per contig reference coverage" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("test_contig_nc1.fa")
          target = Tester.testpath("test_reference_nc1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          comp.per_query_contig_reference_coverage
          assert_equal (2/3.0), comp.assembly["contig1"].reference_coverage
        end
      end
    end

    should "01-1a get per contig reference coverage on protein" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("test_contig_nc1.fa")
          target = Tester.testpath("test_reference_aa1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          comp.per_query_contig_reference_coverage
          assert_equal (2/3.0), comp.assembly["contig1"].reference_coverage
        end
      end
    end

    should "01e raise error because you can't have protein queries" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("test_reference_aa1.fa")
          target = Tester.testpath("test_contig_nc1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          assert_raises Transrate::TransrateError do
            comp.per_query_contig_reference_coverage
          end
        end
      end
    end

    should "02-2n calculate coverage for each reference sequence" do
      Dir.mktmpdir do |tmpdir|
      # tmpdir = Dir.mktmpdir
      # puts tmpdir
        Dir.chdir tmpdir do
          query = Tester.testpath("test_contig_nc2.fa")
          target = Tester.testpath("test_reference_nc1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          comp.per_query_contig_reference_coverage
          comp.per_target_contig_reference_coverage crbblast
          # answer should be 290/300.0
          assert_equal 29/30.0, comp.reference["reference1"].reference_coverage
        end
      end
    end

    should "02-3n calculate coverage for each reference sequence" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("test_contig_nc3.fa")
          target = Tester.testpath("test_reference_nc1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          crbblast.add_missing
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          comp.per_query_contig_reference_coverage
          comp.per_target_contig_reference_coverage crbblast
          # answer should be 1.0000
          assert_equal 1.00, comp.reference["reference1"].reference_coverage
        end
      end
    end

    should "02-3a calculate coverage for each reference sequence" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("test_contig_nc3.fa")
          target = Tester.testpath("test_reference_aa1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          crbblast.add_missing
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          comp.per_query_contig_reference_coverage
          comp.per_target_contig_reference_coverage crbblast
          # answer should be 1.0000
          assert_equal 1.00, comp.reference["reference2"].reference_coverage
        end
      end
    end

    should "02-4n calculate coverage for each reference sequence" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("test_contig_nc4.fa")
          target = Tester.testpath("test_reference_nc1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          comp.per_query_contig_reference_coverage
          comp.per_target_contig_reference_coverage crbblast
          assert_equal 0.5, comp.reference["reference1"].reference_coverage
        end
      end
    end

    should "02-4a calculate coverage for each reference sequence" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("test_contig_nc4.fa")
          target = Tester.testpath("test_reference_aa1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          comp.per_query_contig_reference_coverage
          comp.per_target_contig_reference_coverage crbblast
          assert_equal 0.5, comp.reference["reference2"].reference_coverage
        end
      end
    end

    should "02-5a calculate coverage for each reference sequence" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("test_contig_nc5.fa")
          target = Tester.testpath("test_reference_aa1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          crbblast.add_missing
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          comp.per_query_contig_reference_coverage
          comp.per_target_contig_reference_coverage crbblast
          assert_equal (2/3.0), comp.reference["reference2"].reference_coverage
        end
      end
    end

    should "02-5n calculate coverage for each reference sequence" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("test_contig_nc5.fa")
          target = Tester.testpath("test_reference_nc1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          crbblast.add_missing
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          comp.per_query_contig_reference_coverage
          comp.per_target_contig_reference_coverage crbblast
          assert_equal (2/3.0), comp.reference["reference1"].reference_coverage
        end
      end
    end

    should "02-6a calculate coverage for each reference sequence" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("test_contig_nc6.fa")
          target = Tester.testpath("test_reference_aa1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          crbblast.add_missing
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          comp.per_query_contig_reference_coverage
          comp.per_target_contig_reference_coverage crbblast
          assert_equal (1/3.0), comp.reference["reference2"].reference_coverage
        end
      end
    end

    should "02-6n calculate coverage for each reference sequence" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("test_contig_nc6.fa")
          target = Tester.testpath("test_reference_nc1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          crbblast.add_missing
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          comp.per_query_contig_reference_coverage
          comp.per_target_contig_reference_coverage crbblast
          assert_equal (1/3.0), comp.reference["reference1"].reference_coverage
        end
      end
    end

    should "02-7a calculate coverage for each reference sequence" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("test_contig_nc7.fa")
          target = Tester.testpath("test_reference_aa1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          crbblast.add_missing
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          comp.per_query_contig_reference_coverage
          comp.per_target_contig_reference_coverage crbblast
          assert_equal (1/3.0), comp.reference["reference2"].reference_coverage
        end
      end
    end

    should "02-7n calculate coverage for each reference sequence" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("test_contig_nc7.fa")
          target = Tester.testpath("test_reference_nc1.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          crbblast.add_missing
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          comp.per_query_contig_reference_coverage
          comp.per_target_contig_reference_coverage crbblast
          assert_equal (1/3.0), comp.reference["reference1"].reference_coverage
        end
      end
    end

    should "03 calculate all metrics" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          query = Tester.testpath("assembly.2.fa")
          target = Tester.testpath("Os.protein.2.fa")
          crbblast = CRB_Blast::CRB_Blast.new query, target
          crbblast.run(1e-5, 1, true)
          assembly = Transrate::Assembly.new(query)
          reference = Transrate::Assembly.new(target)
          comp = Transrate::ComparativeMetrics.new(assembly, reference, 1)
          comp.get_reference_hits crbblast
          comp.per_query_contig_reference_coverage
          comp.per_target_contig_reference_coverage crbblast
          assert_equal 11, comp.comp_stats[:CRBB_hits], "CRBB hits"
          assert_equal 11, comp.comp_stats[:n_contigs_with_CRBB], "n_contigs_with_CRBB"
          assert_equal 0.84615, comp.comp_stats[:p_contigs_with_CRBB].round(5), "p_contigs_with_CRBB"
          assert_equal 0.55, comp.comp_stats[:rbh_per_reference], "rbh_per_reference"
          assert_equal 10, comp.comp_stats[:n_refs_with_CRBB], "n_refs_with_CRBB"
          assert_equal 0.5, comp.comp_stats[:p_refs_with_CRBB], "p_refs_with_CRBB"
          assert_equal 10, comp.comp_stats[:cov25], "cov25"
          assert_equal 10, comp.comp_stats[:cov50], "cov50"
          assert_equal 7, comp.comp_stats[:cov75], "cov75"
          assert_equal 6, comp.comp_stats[:cov85], "cov85"
          assert_equal 3, comp.comp_stats[:cov95], "cov95"
          assert_equal 0.5, comp.comp_stats[:p_cov25], "p_cov25"
          assert_equal 0.5, comp.comp_stats[:p_cov50], "p_cov50"
          assert_equal 0.35, comp.comp_stats[:p_cov75], "p_cov75"
          assert_equal 0.3, comp.comp_stats[:p_cov85], "p_cov85"
          assert_equal 0.15, comp.comp_stats[:p_cov95], "p_cov95"
          assert_equal 0.37261, comp.comp_stats[:reference_coverage].round(5), "reference_coverage"

        end
      end
    end



  end
end
