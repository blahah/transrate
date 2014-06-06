#!/usr/bin/env	ruby

require 'helper'

class TestContigMetrics < Test::Unit::TestCase

  context "transrate" do

    setup do
      query = "test/assembly.fasta"
      # query = "test/clg.assembly.fa"
      assembly = Transrate::Assembly.new(query)
      @contig_metrics = Transrate::ContigMetrics.new(assembly)
    end

    should "run metrics on assembly" do
      @contig_metrics.run
      assert @contig_metrics.has_run
    end

    should "get gc content" do
      @contig_metrics.run
      assert_equal @contig_metrics.gc.round(5), 0.37672
    end

    should "get gc skew" do
      @contig_metrics.run
      assert_equal @contig_metrics.gc_skew.round(5), 0.00440
    end

    should "get at skew" do
      @contig_metrics.run
      assert_equal @contig_metrics.at_skew.round(5), -0.00718
    end

    should "get CpG density" do
      @contig_metrics.run
      assert_equal @contig_metrics.cpg.round(5), 0.52828
    end

    should "get linguistic complexity" do
      @contig_metrics.run
      assert_equal @contig_metrics.linguistic_complexity.round(5), 0.26599
    end

    should "get the number and proportion of Ns" do
      @contig_metrics.run
      assert_equal @contig_metrics.bases_n, 2
      assert_equal @contig_metrics.proportion_n.round(5), 0.00033
    end
  end
end