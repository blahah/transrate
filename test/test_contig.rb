require 'helper'
require 'bio'

class TestContig < Test::Unit::TestCase

  context "Contig" do

    setup do
      seq = Bio::Sequence.new 'ATGCGTGTATATACGCGTAG'
      @contig = Transrate::Contig.new seq
    end

    should "know the number and proportion of each base it contains" do
      assert_equal 5, @contig.bases_a, "count of base a"
      assert_equal 0.25, @contig.prop_a, "proportion of base a"
      assert_equal 3, @contig.bases_c, "count of base c"
      assert_equal 0.15, @contig.prop_c, "proportion of base c"
      assert_equal 6, @contig.bases_g, "count of base g"
      assert_equal 0.3, @contig.prop_g, "proportion of base g"
      assert_equal 6, @contig.bases_t, "count of base t"
      assert_equal 0.3, @contig.prop_t, "proportion of base t"
      assert_equal 0, @contig.bases_n, "count of base n"
      assert_equal 0.0, @contig.prop_n, "proportion of base n"
    end

    should "know how many of each two-base pair it contains" do
      assert_equal 3, @contig.dibase_composition[:cg], "cg count"
      assert_equal 3, @contig.dibase_composition[:at], "at count"
      assert_equal 2, @contig.dibase_composition[:tg], "tg count"
    end

    should "know its own gc content" do
      assert_equal 9, @contig.bases_gc, "count of bases that are c or g"
      assert_equal 0.45, @contig.prop_gc.round(2),
                   "proportion of bases that are c or g"
    end

    should "know its own base-pair skew" do
      assert_equal 0.33, @contig.gc_skew.round(2), "gc skew"
      assert_equal -0.09, @contig.at_skew.round(2), "at skew"
    end

    should "know its own CpG count and density" do
      assert_equal 3, @contig.cpg_count, "cpg count"
      assert_equal 66.67, @contig.cpg_ratio.round(2), "cpg ratio"
    end

    should "know the length of its own longest orf" do
      assert_equal 6, @contig.orf_length, "orf length"
    end


    should "calculate linguistic complexity for a long sequence" do
      alphabet = ["A", "C", "G", "T"]
      seq = ""
      50000.times do
        seq << alphabet.sample
      end
      seq = Bio::Sequence.new seq
      contig = Transrate::Contig.new seq
      assert contig.linguistic_complexity(6) > 0.98, "linguistic complexity"
    end

    should "know its own linguistic complexity" do
      assert_equal 0.0586, @contig.linguistic_complexity(4).round(4),
                   "linguistic complexity k=4"
      assert_equal 0.0037, @contig.linguistic_complexity(6).round(4),
                   "linguistic complexity k=6"
    end

  end
end
