require 'helper'
require 'tmpdir'

class TestReadMetrics < Test::Unit::TestCase

  context "ReadMetrics" do

    setup do
      query = File.join(File.dirname(__FILE__), 'data',
                        'sorghum_transcript.fa')
      assembly = Transrate::Assembly.new(query)
      @read_metrics = Transrate::ReadMetrics.new(assembly)
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
          assert_equal 22.91, stats[:mean_coverage].round(2), 'mean coverage'
          assert_equal 11, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.007,
                       stats[:p_uncovered_bases].round(3),
                       'p uncovered bases'
        end
      end
    end

    should "find read pairs that support scaffolding" do
      left = File.join(File.dirname(__FILE__), 'data', 'bridging_reads.l.fastq')
      right = File.join(File.dirname(__FILE__),
                        'data',
                        'bridging_reads.r.fastq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left, right)
          stats = @read_metrics.read_stats
          assert_equal 1, stats[:potential_bridges], 'potential bridges'
        end
      end
    end

    should "count per-base coverage" do

    end

    should "find median coverage" do
    end

    should "identify contigs with at least one uncovered base" do
    end

    should "identify contigs with coverage <1" do
    end

    should "identify contigs with coverage <10" do
    end

  end

end
