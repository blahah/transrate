require 'helper'

class TestContigMetrics < Test::Unit::TestCase

  context "transrate" do

    setup do
      querypath = File.join(File.dirname(__FILE__), 'data',
                            'assembly.fasta')
      assembly = Transrate::Assembly.new(querypath)
      @contig_metrics = Transrate::ContigMetrics.new(assembly)
    end

    should "run metrics on assembly" do
      @contig_metrics.run
      assert @contig_metrics.has_run
    end

    should "get gc content" do
      @contig_metrics.run
      assert_equal 0.37672, @contig_metrics.gc_prop.round(5)
    end

    should "get gc skew" do
      @contig_metrics.run
      assert_equal 0.00440, @contig_metrics.gc_skew.round(5)
    end

    should "get at skew" do
      @contig_metrics.run
      assert_equal -0.00718, @contig_metrics.at_skew.round(5)
    end

    should "get CpG density" do
      @contig_metrics.run
      assert_equal 1.51939, @contig_metrics.cpg_ratio.round(5)
    end

    should "get linguistic complexity" do
      @contig_metrics.run
      assert_equal 0.26526, @contig_metrics.linguistic_complexity.round(5)
    end

    should "get the number and proportion of Ns" do
      @contig_metrics.run
      assert_equal 2, @contig_metrics.bases_n
      assert_equal 0.00033, @contig_metrics.proportion_n.round(5)
    end
  end
end
