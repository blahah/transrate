require 'helper'
require 'bio'
require 'benchmark'

class TestContig < MiniTest::Test

  context "Contig" do

    setup do
      seq = Bio::FastaFormat.new ">test;\nATGCGTGTATATACGCGTAG" # cg=3, gc=2, c*g=20
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

    should "know the length of its own longest orf" do
      assert_equal 6, @contig.orf_length, "orf length"
    end

    should "not break when there is a null byte in the sequence" do
      seq = Bio::FastaFormat.new ">test\nAAAAAAAAAAAA\0AAAAAAAAAAA"
      contig = Transrate::Contig.new seq
      assert_equal 7, contig.orf_length, "orf length"
    end

    should "get orf length from ATG to stop" do
      s = "AAATAGATAGATAGAATGAATCAACTAACTAAATAA"
      # seq = Bio::FastaFormat.new ">test\n#{s}\n>test_r\n#{s.reverse}\n"
      seq = Bio::FastaFormat.new ">test\n#{s}\n>test_r"
      contig = Transrate::Contig.new seq
      assert_equal 6, contig.orf_length, "orf length"
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

    should "classify contig" do
      assert_equal :bad, @contig.classify(0.5), "contig is not bad"
    end

    should "strip trailing semicolons from FASTA entry IDs" do
      assert_equal "test", @contig.name, "trailing semicolon was not stripped"
    end

  end
end
