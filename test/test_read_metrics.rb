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

    should "fail with paired reads that are non-sequential" do
      test={}
      testsam = File.join(File.dirname(__FILE__), 'data', 'nonsequential.sam')
      assert_raise RuntimeError do
        @read_metrics.analyse_read_mappings(testsam, 200, 50, true)
      end
    end

    should "calculate read mapping statistics with paired and unpaired input" do
      left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      unpaired = File.join(File.dirname(__FILE__), 'data', '150uncovered.u.fq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left, right, unpaired, nil)
          #@read_metrics.run(left)
          stats = @read_metrics.read_stats
          assert @read_metrics.has_run
          # checking sam file is consistent
          total=0
          bases=0
          allbases = 0
          samfile = "150uncovered.l.fq.150uncovered.r.fq.150uncovered.u.fq.sorghum_transcript.sam"
          File.open(samfile).each_line do |line|
            unless line[0] == "@"
              allbases += line.split("\t")[9].size
            end
            if line=~/NM:i:([0-9]+)/
              total += $1.to_i
              seq = line.split("\t")[9]
              bases += seq.length
            end
          end
          assert_equal 39996, bases, "total number of bases from sam"
          assert_equal 399, total, "sum of NM from sam"
          #
          assert_equal allbases, @read_metrics.total_bases, "total bases are equal"
          assert_equal 42152, @read_metrics.total_bases, "total bases"
          assert_equal 223, stats[:num_pairs], 'number of read pairs'
          assert_equal 222, stats[:mapped_pairs], 'number mapping'
          assert_equal 94.90, stats[:percent_mapping].round(2),
                       'percent mapping'
          assert_equal 215, stats[:good_mappings], 'good mapping'
          assert_equal 96.41, stats[:pc_good_mapping].round(2),'percent good mapping'
          assert_equal 7, stats[:bad_mappings], 'bad mapping'
          assert_equal 24.2, stats[:mean_coverage].round(2), 'mean coverage'
          assert_equal 12, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.008,
                       stats[:p_uncovered_bases].round(3),
                       'p uncovered bases'
          assert_equal 0.00947, stats[:edit_distance_per_base].round(5),
                       'edit distance'
          assert_equal 63.66636,stats[:coverage_variance].round(5),
                       'coverage variance'
        end
      end
    end
    should "calculate read mapping statistics with paired input" do
      left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left, right, nil, nil)
          stats = @read_metrics.read_stats
          assert @read_metrics.has_run
          # checking sam file is consistent
          total=0
          bases=0
          allbases = 0
          samfile = "150uncovered.l.fq.150uncovered.r.fq..sorghum_transcript.sam"
          File.open(samfile).each_line do |line|
            unless line[0] == "@"
              allbases += line.split("\t")[9].size
            end
            if line=~/NM:i:([0-9]+)/
              total += $1.to_i
              seq = line.split("\t")[9]
              bases += seq.length
            end
          end
          assert_equal 37921, bases, "total number of bases from sam"
          assert_equal 387, total, "sum of NM from sam"
          #
          assert_equal allbases, @read_metrics.total_bases, "total bases are equal"
          assert_equal 39909, @read_metrics.total_bases, "total bases"
          assert_equal 223, stats[:num_pairs], 'number of read pairs'
          assert_equal 222, stats[:mapped_pairs], 'number mapping'
          assert_equal 95.07, stats[:percent_mapping].round(2),
                       'percent mapping'
          assert_equal 215, stats[:good_mappings], 'good mapping'
          assert_equal 96.41,
                       stats[:pc_good_mapping].round(2),
                       'percent good mapping'
          assert_equal 7, stats[:bad_mappings], 'bad mapping'
          assert_equal 22.86, stats[:mean_coverage].round(2), 'mean coverage'
          assert_equal 12, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.008,
                       stats[:p_uncovered_bases].round(3),
                       'p uncovered bases'
          assert_equal 0.0097, stats[:edit_distance_per_base].round(5),
                       'edit distance'
          assert_equal 57.94268,stats[:coverage_variance].round(5),
                       'coverage variance'
        end
      end
    end
    should "calculate read mapping statistics for unpaired input" do
      unpaired = File.join(File.dirname(__FILE__), 'data', '150uncovered.u.fq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(nil,nil,unpaired,nil)
          stats = @read_metrics.read_stats
          assert @read_metrics.has_run
          # checking sam file is consistent
          total=0
          bases=0
          allbases=0
          samfile = "..150uncovered.u.fq.sorghum_transcript.sam"
          File.open(samfile).each_line do |line|
            unless line[0] == "@"
              allbases += line.split("\t")[9].size
            end
            if line=~/NM:i:([0-9]+)/
              total += $1.to_i
              seq = line.split("\t")[9]
              bases += seq.length
            end
          end
          assert_equal 2075, bases, "total number of bases from sam"
          assert_equal 12, total, "sum of NM from sam"
          #
          assert_equal allbases, @read_metrics.total_bases, "total bases are equal"
          assert_equal 2243, @read_metrics.total_bases, "total bases"
          assert_equal 0, stats[:num_pairs], 'number of read pairs'
          assert_equal 23, stats[:mapped_single_reads], 'number mapping'
          assert_equal 92.0, stats[:percent_mapping].round(2),
                       'percent mapping'
          assert_equal 0, stats[:good_mappings], 'good mapping'
          assert_equal 0, stats[:bad_mappings], 'bad mapping'
          assert_equal 1.33, stats[:mean_coverage].round(2), 'mean coverage'
          assert_equal 691, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.444,
                       stats[:p_uncovered_bases].round(3),
                       'p uncovered bases'
          assert_equal 0.00535, stats[:edit_distance_per_base].round(5),
                       'edit distance'
          assert_equal 1.12367,stats[:coverage_variance].round(5),
                       'coverage variance'
        end
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
          assert_equal 10.84, stats[:mean_coverage].round(2), 'mean coverage'
          assert_equal 53, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.034,
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
          assert_equal 919, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.591,
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
          assert_equal 919, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.591,
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
          assert_equal 919, stats[:n_uncovered_bases], 'n uncovered bases'
          assert_equal 0.591,
                       stats[:p_uncovered_bases].round(3),
                       'p uncovered bases'
        end
      end
    end

    should "get metrics for individual contigs" do
      left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @read_metrics.run(left, right, nil, nil)
          contigs = []
          @assembly.each do |name, contig|
            contigs << contig
          end
          a = contigs[0].read_metrics
          b = contigs[1].read_metrics
          edit_a = a[:edit_distance_per_base].round(5)
          edit_b = b[:edit_distance_per_base].round(5)
          assert_equal 0.00801, edit_a, "edit distance"
          assert_equal 0.01081, edit_b, "edit distance"
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
          @read_metrics.run(left, right, nil, nil)
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
          @read_metrics.run(left.join(","), right.join(","), nil, nil)
          stats = @read_metrics.read_stats
          assert_equal 228, stats[:num_pairs], 'number of read pairs'
        end
      end
    end

  end
end
