#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestTransrater < Test::Unit::TestCase

  context "transrater" do

    setup do
      @assembly = File.join(File.dirname(__FILE__), 'data', 'assembly.2.fa')
      @reference = File.join(File.dirname(__FILE__), 'data', 'Os.protein.2.fa')
      @test_assembly = File.join(File.dirname(__FILE__), 'data', 'test_ref_align_assembly.fa')
      @test_genome = File.join(File.dirname(__FILE__), 'data', 'test_ref_align_genome.fa')
      @rater = Transrate::Transrater.new(@assembly, @reference, nil)
      @rater_aligner = Transrate::Transrater.new(@test_assembly, nil, @test_genome)
      @left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      @right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
    end

    should "create transrater object" do
      assert @rater
    end

    should "raise error when assembly input is nil" do
      assert_raise RuntimeError do
        rater = Transrate::Transrater.new(nil, @reference, nil)
      end
    end

    should "handle assembly as an assemby object" do
      assembly_object = Transrate::Assembly.new(@assembly)
      rater = Transrate::Transrater.new(assembly_object, @reference, nil)
    end

    should "handle reference as an assemby object" do
      reference_object = Transrate::Assembly.new(@reference)
      rater = Transrate::Transrater.new(@assembly, reference_object, nil)
    end

    should "handle genome as an assembly object" do
      genome_object = Transrate::Assembly.new(@test_genome)
      rater = Transrate::Transrater.new(@test_assembly, nil, genome_object)
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

    should "run read metrics with input fastq files" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          stats = @rater.read_metrics(@left, @right)
          assert_equal 223, stats.read_stats[:num_pairs]
        end
      end
    end

    should "get assembly score" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          all = @rater.all_metrics(@left, @right)
          score = @rater.assembly_score
          assert_equal 0.23282, score.round(5)
        end
      end
    end

  end
end
