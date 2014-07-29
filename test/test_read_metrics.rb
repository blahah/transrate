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
          assert @read_metrics.has_run
          assert_equal 223, stats[:num_pairs], 'number of read pairs'
          assert_equal 202, stats[:total_mappings], 'number mapping'
          assert_equal 90.58, stats[:percent_mapping].round(2),
                       'percent mapping'
          assert_equal 202, stats[:good_mappings], 'good mapping'
          assert_equal 90.58,
                       stats[:pc_good_mapping].round(2),
                       'percent good mapping'
          assert_equal 0, stats[:bad_mappings], 'bad mapping'
          assert_equal 22.86, stats[:mean_coverage].round(2), 'mean coverage'
          assert_equal 12, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.008,
                       stats[:p_uncovered_bases].round(3),
                       'p uncovered bases'
          assert_equal 0.01009, stats[:edit_distance_per_base].round(5),
                       'edit distance'
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
          a = contigs[0].read_metrics[:edit_distance_per_base].round(5)
          b = contigs[1].read_metrics[:edit_distance_per_base].round(5)
          assert_equal 0.00788, a, "edit distance"
          assert_equal 0.01157, b, "edit distance"
          a = contigs[0].read_metrics[:low_uniqueness_bases]
          b = contigs[1].read_metrics[:low_uniqueness_bases]
          assert_equal 0, a, "low uniqueness bases"
          assert_equal 0, b, "low uniqueness bases"
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

  end

end
