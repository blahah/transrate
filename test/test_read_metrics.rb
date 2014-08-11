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
          # checking sam file is consistent
          total=0
          bases=0
          File.open("150uncovered.l.fq.150uncovered.r.fq.sorghum_transcript.sam").each_line do |line|
            if line=~/NM:i:([0-9]+)/
              total += $1.to_i
              seq = line.split("\t")[9]
              bases += seq.length
            end

          end
          assert_equal 37921, bases, "total number of bases from sam"
          assert_equal 387, total, "sum of NM from sam"
          #
          assert_equal 37921, @read_metrics.total_bases, "total bases"
          assert_equal 223,     stats[:num_pairs], 'number of read pairs'
          assert_equal 202,     stats[:total_mappings], 'number mapping'
          assert_equal 90.58,   stats[:percent_mapping].round(2),
                       'percent mapping'
          assert_equal 202,     stats[:good_mappings], 'good mapping'
          assert_equal 90.58,   stats[:pc_good_mapping].round(2),
                       'percent good mapping'
          assert_equal 0,       stats[:bad_mappings], 'bad mapping'
          assert_equal 22.86,   stats[:mean_coverage].round(2), 'mean coverage'
          assert_equal 12,      stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.008,   stats[:p_uncovered_bases].round(3),
                       'p uncovered bases'

          assert_equal 0.01021, stats[:edit_distance_per_base].round(5),
                       'edit distance'
          assert_equal 57.94268,stats[:coverage_variance].round(5),
                       'coverage variance'
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
          edit_a = a[:edit_distance_per_base].round(5)
          edit_b = b[:edit_distance_per_base].round(5)
          assert_equal 0.00832, edit_a, "edit distance"
          assert_equal 0.01138, edit_b, "edit distance"
          uniq_a = a[:low_uniqueness_bases]
          uniq_b = a[:low_uniqueness_bases]
          assert_equal 0, uniq_a, "low uniqueness bases"
          assert_equal 0, uniq_b, "low uniqueness bases"
          var_a = contigs[0].effective_variance
          var_b = contigs[1].effective_variance
          assert_equal 109.75611, var_a.round(5)
          assert_equal 13.69750, var_b.round(5)
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

  end

end
