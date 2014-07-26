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
          assert_equal 222, stats[:mapped_pairs], 'number mapping'
          assert_equal 94.90, stats[:percent_mapping].round(2),
                       'percent mapping'
          assert_equal 215, stats[:good_mappings], 'good mapping'
          assert_equal 96.41, stats[:pc_good_mapping].round(2),'percent good mapping'
          assert_equal 7, stats[:bad_mappings], 'bad mapping'
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
          @read_metrics.run(left, right, nil)
          stats = @read_metrics.read_stats
          assert @read_metrics.has_run
          assert_equal 223, stats[:num_pairs], 'number of read pairs'
          assert_equal 222, stats[:mapped_pairs], 'number mapping'
          assert_equal 95.07, stats[:percent_mapping].round(2),
                       'percent mapping'
          assert_equal 215, stats[:good_mappings], 'good mapping'
          assert_equal 96.41,
                       stats[:pc_good_mapping].round(2),
                       'percent good mapping'
          assert_equal 7, stats[:bad_mappings], 'bad mapping'
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
          assert_equal 23, stats[:mapped_single_reads], 'number mapping'
          assert_equal 92.0, stats[:percent_mapping].round(2),
                       'percent mapping'
          assert_equal 0, stats[:good_mappings], 'good mapping'
          assert_equal 0, stats[:bad_mappings], 'bad mapping'
          assert_equal 1.34, stats[:mean_coverage].round(2), 'mean coverage'
          assert_equal 689, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.443,
                       stats[:p_uncovered_bases].round(3),
                       'p uncovered bases'
        end
      end
    end

    should "fail with paired reads that are non-sequential" do
      test={}
      testsam = File.join(File.dirname(__FILE__), 'data', 'nonsequential.sam')
      assert_raise RuntimeError do
        @read_metrics.analyse_read_mappings(testsam, 200, 50, true)
      end
    end

    should "run with fr strand specific input" do
      left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      unpaired = File.join(File.dirname(__FILE__), 'data', '150uncovered.u.fq')
      fr = "fr"
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left,right,unpaired,fr)
          stats = @read_metrics.read_stats
          assert @read_metrics.has_run
          assert_equal 223, stats[:num_pairs], 'number of read pairs'
          assert_equal 103, stats[:mapped_pairs], 'number mapping'
          assert_equal 43.1, stats[:percent_mapping].round(2),
                       'percent mapping'
          assert_equal 99, stats[:good_mappings], 'good mapping'
          assert_equal 44.39, stats[:pc_good_mapping].round(2),'percent good mapping'
          assert_equal 4, stats[:bad_mappings], 'bad mapping'
          assert_equal 10.86, stats[:mean_coverage].round(2), 'mean coverage'
          assert_equal 51, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.033,
                       stats[:p_uncovered_bases].round(3),
                       'p uncovered bases'
        end
      end
    end

    should "run with rf strand specific input" do
      left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      unpaired = File.join(File.dirname(__FILE__), 'data', '150uncovered.u.fq')
      rf = "rf"
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left,right,unpaired,rf)
          stats = @read_metrics.read_stats
          assert @read_metrics.has_run
          assert_equal 471, stats[:num_reads], 'number of reads'
          assert_equal 223, stats[:num_pairs], 'number of read pairs'
          assert_equal 119, stats[:mapped_pairs], 'number mapping'
          assert_equal 51.17, stats[:percent_mapping].round(2),
                       'percent mapping'
          assert_equal 67, stats[:good_mappings], 'good mapping'
          assert_equal 30.04, stats[:pc_good_mapping].round(2),'percent good mapping'
          assert_equal 52, stats[:bad_mappings], 'bad mapping'
          assert_equal 0.57, stats[:mean_coverage].round(2), 'mean coverage'
          assert_equal 917, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.59,
                       stats[:p_uncovered_bases].round(3),
                       'p uncovered bases'
        end
      end
    end

    should "run with ff strand specific input" do
      left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      unpaired = File.join(File.dirname(__FILE__), 'data', '150uncovered.u.fq')
      rf = "ff"
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left,right,unpaired,rf)
          stats = @read_metrics.read_stats
          assert @read_metrics.has_run
          assert_equal 471, stats[:num_reads], 'number of reads'
          assert_equal 223, stats[:num_pairs], 'number of read pairs'
          assert_equal 209, stats[:mapped_pairs], 'number mapping'
          assert_equal 46.5, stats[:percent_mapping].round(2),
                       'percent mapping'
          assert_equal 0, stats[:good_mappings], 'good mapping'
          assert_equal 209, stats[:bad_mappings], 'bad mapping'
          assert_equal 0.57, stats[:mean_coverage].round(2), 'mean coverage'
          assert_equal 917, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.59,
                       stats[:p_uncovered_bases].round(3),
                       'p uncovered bases'
        end
      end
    end

    should "run with f strand specific input" do
      unpaired = File.join(File.dirname(__FILE__), 'data', '150uncovered.u.fq')
      f = "f"
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(nil,nil,unpaired,f)
          stats = @read_metrics.read_stats
          assert @read_metrics.has_run
          assert_equal 25, stats[:num_reads], 'number of reads'
          assert_equal 0, stats[:mapped_pairs], 'number of mapping pairs'
          assert_equal 10, stats[:total_mapped_reads], "mapped reads"
          assert_equal 40.0, stats[:percent_mapping].round(2),
                       'percent mapping'
          assert_equal 0, stats[:good_mappings], 'good mapping'
          assert_equal 0, stats[:bad_mappings], 'bad mapping'
          assert_equal 0.57, stats[:mean_coverage].round(2), 'mean coverage'
          assert_equal 917, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.59,
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
          @read_metrics.run(left, right, nil)
          stats = @read_metrics.read_stats
          assert_equal 1, stats[:potential_bridges], 'potential bridges'
        end
      end
    end
  end
end
