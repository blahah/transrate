#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestTransrater < MiniTest::Test

  context "transrater" do

    setup do
      @assembly = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.fa')
      @reference = File.join(File.dirname(__FILE__), 'data', 'Os.protein.2.fa')
      @rater = Transrate::Transrater.new(@assembly, @reference)
      @left = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.1.fastq')
      @right = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.2.fastq')
    end

    should "create transrater object" do
      assert @rater
    end

    should "raise error when assembly input is nil" do
      assert_raises Transrate::TransrateError do
        rater = Transrate::Transrater.new(nil, @reference)
      end
    end

    should "handle assembly as an assembly object" do
      assembly_object = Transrate::Assembly.new(@assembly)
      rater = Transrate::Transrater.new(assembly_object, @reference)
    end

    should "handle reference as an assembly object" do
      reference_object = Transrate::Assembly.new(@reference)
      rater = Transrate::Transrater.new(@assembly, reference_object)
    end

    should "run assembly metrics without input fastq files" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          stats = @rater.assembly_metrics
          assert_equal 1692, stats.n50
          assert_equal 137748, stats.n_bases
        end
      end
    end

    should "run read metrics with input fastq files" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          stats = @rater.read_metrics(@left, @right)
          assert_equal 25006, stats.read_stats[:fragments]
        end
      end
    end

    should "get assembly score" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          all = @rater.all_metrics(@left, @right)
          score = @rater.assembly_score
          assert_equal 0.154, score.round(3) # regression test
        end
      end
    end

  end
end
