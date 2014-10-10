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
          assert_equal 39196,   @read_metrics.total_bases, "total bases"
          assert_equal 223,     stats[:num_pairs], 'number of read pairs'
          assert_equal 215,     stats[:total_mappings], 'number mapping'
          assert_equal 96.41,   stats[:percent_mapping].round(2),
                       'percent mapping'
          assert_equal 144,     stats[:good_mappings], 'good mapping'
          assert_equal 64.57,   stats[:pc_good_mapping].round(2),
                       'percent good mapping'
          assert_equal 79,      stats[:bad_mappings], 'bad mapping'
          assert_equal 52.84,   stats[:mean_mapq].round(2), 'mean mapq'
          assert_equal 2, stats[:potential_bridges], 'bridges'
          assert_equal 2, stats[:n_uncovered_base_contigs], 'uncovered base contig'
          assert_equal 0, stats[:n_uncovered_contigs], 'uncovered contig'
          assert_equal 0, stats[:n_lowcovered_contigs], 'lowcovered contig'
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

          edit_a = a[:inverse_edit_dist].round(5)
          edit_b = b[:inverse_edit_dist].round(5)
          assert_equal 0.99507, edit_a, "edit distance 1"
          assert_equal 0.98987, edit_b, "edit distance 2"

          uniq_a = a[:p_unique_bases]
          uniq_b = a[:p_unique_bases]
          assert_equal 1, uniq_a, "unique bases"
          assert_equal 1, uniq_b, "unique bases"

          assert_equal 0.68539, a[:p_good].round(5), "proportion of good mappings"
          assert_equal 0.62406, b[:p_good].round(5), "proportion of good mappings"

          # uncovered bases
          unc_a = contigs[0].uncovered_bases
          unc_b = contigs[1].uncovered_bases
          assert_equal 11, unc_a, "uncovered bases"
          assert_equal 1, unc_b, "uncovered bases"

          prop_unc_a = a[:p_bases_covered]
          prop_unc_b = b[:p_bases_covered]
          assert_equal 0.98497, prop_unc_a.round(5), "proportion covered bases"
          assert_equal 0.99878, prop_unc_b.round(5), "proportion covered bases"

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
          assert_equal 0.28917, edit_a, "probability not segmented 1"
          assert_equal 0.79553, edit_b, "probability not segmented 2"

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
          assert_equal 228, stats[:num_pairs], 'number of read pairs'
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
