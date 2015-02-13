require 'helper'
require 'tmpdir'

class TestReadMetrics < Test::Unit::TestCase

  context "ReadMetrics" do

    setup do
      query = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.fa')
      @assembly = Transrate::Assembly.new(query)
      @read_metrics = Transrate::ReadMetrics.new(@assembly)
    end

    teardown do
      if File.exist?("test/data/sorghum_100.fa.fai")
        rm = "rm test/data/sorghum_100.fa.fai"
        `#{rm}`
      end
    end

    should "setup correctly" do
      assert @read_metrics
    end

    should "calculate read mapping statistics" do
      left = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.1.fastq')
      right = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.2.fastq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left, right)
          stats = @read_metrics.read_stats
          assert @read_metrics.has_run, "has run"
          assert_equal 25006,    stats[:fragments], 'number of read pairs'
          assert_equal 21000,    stats[:fragments_mapped], 'number mapping'
          assert_in_delta 0.84,  stats[:p_fragments_mapped].round(4), 0.05
                       'proportion mapping'
          assert_in_delta 17382, stats[:good_mappings], 50, 'good mapping'
          assert_in_delta 0.6945,stats[:p_good_mapping].round(3), 0.005,
                       'percent good mapping'
          assert_in_delta 3637,  stats[:bad_mappings], 50, 'bad mapping'
          assert_equal 1,  stats[:potential_bridges], 'bridges'
          assert_equal 93, stats[:contigs_uncovbase], 'uncovered base contig'
          assert_equal 28, stats[:contigs_uncovered], 'uncovered contig'
          assert_equal 72, stats[:contigs_lowcovered], 'lowcovered contig'
        end
      end
    end

    should "get metrics for individual contigs" do
      left = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.1.fastq')
      right = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.2.fastq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left, right)
          contigs = []
          @assembly.each do |name, contig|
            contigs << contig
          end
          a = contigs[0].read_metrics
          b = contigs[1].read_metrics

          edit_a = a[:p_seq_true].round(5)
          edit_b = b[:p_seq_true].round(5)
          assert_in_delta 0.6597, edit_a.round(4), 0.01, "edit distance 1"
          assert_in_delta 0.9364, edit_b.round(4), 0.01, "edit distance 2"

          assert_in_delta 0.5714, a[:p_good].round(4), 0.01, "prop of good mappings"
          assert_in_delta 0.8302, b[:p_good].round(4), 0.01, "prop of good mappings"

          # uncovered bases
          unc_a = contigs[0].uncovered_bases
          unc_b = contigs[1].uncovered_bases
          assert_in_delta 366, unc_a, 20, "uncovered bases"
          assert_in_delta 64,  unc_b, 20, "uncovered bases"

          prop_unc_a = a[:p_bases_covered]
          prop_unc_b = b[:p_bases_covered]
          assert_in_delta 0.5622, prop_unc_a.round(4), 0.05, "proportion covered bases"
          assert_in_delta 0.9964, prop_unc_b.round(4), 0.05, "proportion covered bases"

        end
      end
    end

    should "classify individual contigs" do
      left = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.1.fastq')
      right = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.2.fastq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left, right)
          contigs = []
          @assembly.each do |name, contig|
            contigs << contig
          end
          assert_equal :bad,  contigs[0].classify(0.5)
          assert_equal :good, contigs[1].classify(0.5)
          assert_equal :bad,  contigs[2].classify(0.5)
          assert_equal :good, contigs[3].classify(0.5)
        end
      end

    end

    should "get segmentation probability for individual contigs" do
      left = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.1.fastq')
      right = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.2.fastq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left, right)
          contigs = []
          @assembly.each do |name, contig|
            contigs << contig
          end
          a = contigs[0].read_metrics
          b = contigs[1].read_metrics

          edit_a = a[:p_not_segmented].round(5)
          edit_b = b[:p_not_segmented].round(5)
          assert_in_delta 0.9965, edit_a.round(4), 0.05, "probability not segmented 1"
          assert_in_delta 0.9644, edit_b.round(4), 0.05, "probability not segmented 2"

        end
      end
    end

    should "find read pairs that support scaffolding" do
      left = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.1.fastq')
      right = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.2.fastq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left, right)
          contigs = []
          @assembly.each do |name, contig|
            if contig.in_bridges > 0
              contigs << contig
            end
          end
          stats = @read_metrics.read_stats
          assert @read_metrics.has_run, "has run"
          assert_equal 1,  stats[:potential_bridges], 'bridges'
          assert_equal "Sb01g002430.1", contigs[0].name
          assert_equal "Sb05g009350.1", contigs[1].name
        end
      end
    end

    should "run on a list of input fastq files" do
      left = []
      right = []
      left << File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right << File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')

      left << File.join(File.dirname(__FILE__), 'data', 'sorghum_100.1.fastq')
      right << File.join(File.dirname(__FILE__), 'data', 'sorghum_100.2.fastq')

      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left.join(","), right.join(","))
          stats = @read_metrics.read_stats
          assert_equal 25229, stats[:fragments], 'number of read pairs'
        end
      end

    end

    should "calculate read length by scanning fastq files" do
      left = []
      right = []
      left << File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right << File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      left << File.join(File.dirname(__FILE__),
                       'data', 'bridging_reads.l.fastq')
      right << File.join(File.dirname(__FILE__),
                        'data', 'bridging_reads.r.fastq')
      @read_metrics.get_read_length(left.join(","), right.join(","))
      assert_equal 100, @read_metrics.read_length
    end

  end

end
