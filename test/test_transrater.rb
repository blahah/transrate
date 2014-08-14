#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestTransrater < Test::Unit::TestCase

  context "transrater" do

    setup do
      @assembly = File.join(File.dirname(__FILE__), 'data', 'assembly.2.fa')
      @reference = File.join(File.dirname(__FILE__), 'data', 'Os.protein.2.fa')
      @rater = Transrate::Transrater.new(@assembly, @reference)
      @left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      @right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      @unpaired = File.join(File.dirname(__FILE__), 'data', '150uncovered.u.fq')
    end

    should "create transrater object" do
      assert @rater
    end

    should "raise error when assembly input is nil" do
      assert_raise RuntimeError do
        rater = Transrate::Transrater.new(nil, @reference)
      end
    end

    should "handle assembly as an assemby object" do
      assembly_object = Transrate::Assembly.new(@assembly)
      rater = Transrate::Transrater.new(assembly_object, @reference)
    end

    should "handle reference as an assemby object" do
      reference_object = Transrate::Assembly.new(@reference)
      rater = Transrate::Transrater.new(@assembly, reference_object)
    end

    should "run assembly metrics without input fastq files" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          stats = @rater.assembly_metrics
          assert_equal 1566, stats.n50
          assert_equal 10331, stats.n_bases
        end
      end
    end

    should "run read metrics with unpaired and paired input" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          stats = @rater.read_metrics(@left, @right, @unpaired)

          assert_equal 223, stats.read_stats[:num_pairs]
          assert_equal 471, stats.read_stats[:num_reads]
        end
      end
    end

    should "run read metrics with paired input" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          stats = @rater.read_metrics(@left,@right)
          assert_equal 223, stats.read_stats[:num_pairs]
          assert_equal 446, stats.read_stats[:num_reads]
        end
      end
    end

    should "run read metrics with unpaired input" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          stats = @rater.read_metrics(nil,nil,@unpaired)
          assert_equal 0, stats.read_stats[:num_pairs]
          assert_equal 25, stats.read_stats[:num_reads]
        end
      end
    end
    should "get assembly score" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          all = @rater.all_metrics(@left, @right)
          score = @rater.assembly_score
          assert_equal 0.18605, score.round(5) # regression test
        end
      end
    end

  end
end
