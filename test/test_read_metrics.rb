require 'helper'
require 'tmpdir'

class TestReadMetrics < Test::Unit::TestCase

  context "ReadMetrics" do

    setup do
      query = File.join(File.dirname(__FILE__), 'data',
                        'sorghum_transcript.fa')
      @assembly = Transrate::Assembly.new(query)
      @read_metrics = Transrate::ReadMetrics.new(@assembly)
    end

    teardown do
      if File.exist?("test/data/sorghum_transcript.fa.fai")
        rm = "rm test/data/sorghum_transcript.fa.fai"
        `#{rm}`
      end
    end

    should "setup correctly" do
      assert @read_metrics
    end

    should "calculate read mapping statistics" do
      left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left, right)
          stats = @read_metrics.read_stats
          assert @read_metrics.has_run, "has run"
          assert_equal 223,     stats[:fragments], 'number of read pairs'
          assert_equal 217,     stats[:fragments_mapped], 'number mapping'
          assert_equal 0.9731,   stats[:p_fragments_mapped].round(4),
                       'proportion mapping'
          assert_equal 129,     stats[:good_mappings], 'good mapping'
          assert_equal 0.5785,   stats[:p_good_mapping].round(4),
                       'percent good mapping'
          assert_equal 88,      stats[:bad_mappings], 'bad mapping'
          assert_equal 2, stats[:potential_bridges], 'bridges'
          assert_equal 2, stats[:contigs_uncovbase], 'uncovered base contig'
          assert_equal 0, stats[:contigs_uncovered], 'uncovered contig'
          assert_equal 0, stats[:contigs_lowcovered], 'lowcovered contig'
          assert_equal 0, stats[:contigs_good], 'good contigs'
        end
      end
    end

    should "get metrics for individual contigs" do
      left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
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
          assert_equal 0.98319, edit_a, "edit distance 1"
          assert_equal 0.96509, edit_b, "edit distance 2"

          assert_equal 0.8046, a[:p_good].round(5),
                       "proportion of good mappings"
          assert_equal 0.875, b[:p_good].round(5), "proportion of good mappings"

          # uncovered bases
          unc_a = contigs[0].uncovered_bases
          unc_b = contigs[1].uncovered_bases
          assert_equal 11, unc_a, "uncovered bases"
          assert_equal 2, unc_b, "uncovered bases"

          prop_unc_a = a[:p_bases_covered]
          prop_unc_b = b[:p_bases_covered]
          assert_equal 0.98497, prop_unc_a.round(5), "proportion covered bases"
          assert_equal 0.99757, prop_unc_b.round(5), "proportion covered bases"

        end
      end
    end

    should "get segmentation probability for individual contigs" do
      left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
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
          assert_equal 0.11649, edit_a, "probability not segmented 1"
          assert_equal 0.85865, edit_b, "probability not segmented 2"

        end
      end
    end

    should "find read pairs that support scaffolding" do
      left = File.join(File.dirname(__FILE__),
                       'data', 'bridging_reads.l.fastq')
      right = File.join(File.dirname(__FILE__),
                        'data', 'bridging_reads.r.fastq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left, right)
          stats = @read_metrics.read_stats
          assert_equal 1, stats[:potential_bridges], 'potential bridges'
        end
      end
    end

    should "run on a list of input fastq files" do
      left = []
      right = []
      left << File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right << File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      left << File.join(File.dirname(__FILE__),
                       'data', 'bridging_reads.l.fastq')
      right << File.join(File.dirname(__FILE__),
                        'data', 'bridging_reads.r.fastq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left.join(","), right.join(","))
          stats = @read_metrics.read_stats
          assert_equal 233, stats[:fragments], 'number of read pairs'
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
