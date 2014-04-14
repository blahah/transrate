#!/usr/bin/env  ruby

require 'helper'

class TestInline < Test::Unit::TestCase

  context "transrate" do

    setup do
      @a = Transrate::Assembly.new("test/assembly.fasta")
      @seq1 = "ATGCCCCTAGGGTAG"
      @seq2 = @seq1.reverse.tr("ACGT", "TGCA")
    end

    should "find longest orf in file" do
      orfs = []
      @a.assembly.each do |entry|
        l = @a.orf_length entry.seq
        orfs << l
      end
      assert_equal orfs.length, 4
      assert_equal orfs, [333, 370, 131, 84]
    end

    should "find longest orf in sequence" do
      l = @a.orf_length(@seq1)
      assert_equal l, 5
    end

  end
end