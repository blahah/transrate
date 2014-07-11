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

    teardown do
      if File.exist?("test/data/sorghum_transcript.fa.fai")
        rm = "rm test/data/sorghum_transcript.fa.fai"
        `#{rm}`
      end
    end

    should "setup correctly" do
      assert @read_metrics
    end

    should "calculate read mapping statistics with paired and unpaired input" do
      left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      unpaired = File.join(File.dirname(__FILE__), 'data', '150uncovered.u.fq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left, right, unpaired)
          #@read_metrics.run(left)
          stats = @read_metrics.read_stats
          assert @read_metrics.has_run
          assert_equal 223, stats[:num_pairs], 'number of read pairs'
          assert_equal 213, stats[:total_mappings], 'number mapping'
          assert_equal 95.52, stats[:percent_mapping].round(2),
                       'percent mapping'
          #assert_equal 205, stats[:good_mappings], 'good mapping'
          #assert_equal 91.93, stats[:pc_good_mapping].round(2),'percent good mapping'
          #assert_equal 1, stats[:bad_mappings], 'bad mapping'
          assert_equal 24.25, stats[:mean_coverage].round(2), 'mean coverage'
          assert_equal 11, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.007,
                       stats[:p_uncovered_bases].round(3),
                       'p uncovered bases'
        end
      end
    end
    should "calculate read mapping statistics with paired input" do
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
    should "calculate read mapping statistics for unpaired input" do
      unpaired = File.join(File.dirname(__FILE__), 'data', '150uncovered.u.fq')

      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(nil,nil,unpaired)
          stats = @read_metrics.read_stats
          assert @read_metrics.has_run
          assert_equal 0, stats[:num_pairs], 'number of read pairs'
          assert_equal 11, stats[:total_mappings], 'number mapping'
          assert_equal Float::INFINITY, stats[:percent_mapping].round(2),
                       'percent mapping'
          #assert_equal 2, stats[:good_mappings], 'good mapping'
          #assert_equal Float::INFINITY,stats[:pc_good_mapping].round(2),'percent good mapping'
          #assert_equal 9, stats[:bad_mappings], 'bad mapping'
          assert_equal 1.34, stats[:mean_coverage].round(2), 'mean coverage'
          assert_equal 689, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.443,
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

  end

end
