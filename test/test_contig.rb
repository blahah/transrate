require 'helper'
require 'bio'
require 'benchmark'

class TestContig < Test::Unit::TestCase

  context "Contig" do

    setup do
      seq = Bio::FastaFormat.new ">test\nATGCGTGTATATACGCGTAG" # cg=3, gc=2, c*g=20
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

    should "calculate dibase composition with ambiguous bases" do
      seq = Bio::FastaFormat.new ">test\nATGGGNCRYTAG"
      contig = Transrate::Contig.new seq
      assert_equal 1, contig.dibase_composition[:at]
      assert_equal 1, contig.dibase_composition[:nn]
      assert_equal 1, contig.dibase_composition[:gn]
    end

    should "benchmark composition" do
      seq = "GCCGTGAGCTTCTTGATCGAGTTCTTCTCCCGCTTCGCGAACGCCTTGGACTCCTNGCACGGG"
      seq << "GTCAGCCCCGCGATGTCGGCGGCCGCGGGCGGGGGG"
      seq = seq * 100000
      seq = Bio::FastaFormat.new ">test\n"+seq
      contig = Transrate::Contig.new seq
      ruby_time = 11 # time taken with the ruby version
      c_time = Benchmark.realtime do |x|
        contig.dibase_composition
      end
      assert c_time*100 < ruby_time, "c faster than ruby"
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

    should "know its own CpG count and ratio" do
      assert_equal 5, @contig.cpg_count, "cpg count"
      assert_equal 5.56, @contig.cpg_ratio.round(2), "cpg ratio"
    end

    should "calculate the CpG ratio" do
      seq = Bio::FastaFormat.new ">test\nAAACGAAA"
      contig = Transrate::Contig.new seq
      assert_equal 8, contig.cpg_ratio, "cpg ratio"
    end

    should "know the length of its own longest orf" do
      assert_equal 6, @contig.orf_length, "orf length"
    end

    should "not break when there is a null byte in the sequence" do
      seq = Bio::FastaFormat.new ">test\nAAAAAAAAAAAA\0AAAAAAAAAAA"
      contig = Transrate::Contig.new seq
      assert_equal 7, contig.orf_length, "orf length"
    end

    should "not fail on bases that aren't ACGTN" do
      seq = Bio::FastaFormat.new ">test\nATGCGTGTARATACGCGTAG"
      contig = Transrate::Contig.new seq
      assert_equal 1, contig.base_composition[:n]
    end

    should "get kmer count with non ACGTN bases" do
      seq = Bio::FastaFormat.new ">test\nATGCGTGTARATACGCGTAG"
      contig = Transrate::Contig.new seq
      assert_equal 0, contig.kmer_count(6, "RRRRRRRRRRRRRRRR")
    end

    should "calculate linguistic complexity for a long sequence" do
      alphabet = ["A", "C", "G", "T"]
      seq = ""
      50000.times do
        seq << alphabet.sample
      end
      seq = Bio::FastaFormat.new ">test\n"+seq
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
