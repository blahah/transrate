#!/usr/bin/env	ruby

require 'helper'

class TestTransrate < Test::Unit::TestCase

  context "transrate" do

    setup do
      @a = Transrate::Assembly.new("test/assembly.fasta")
      @seq1 = "ATGCCCGGGTAG"
    end

    should "run metrics on assembly" do
      ans = @a.run(2) # using 2 threads
      assert_equal ans, true, "should run but returned #{ans}"
    end

    should "find longest orf" do
      len = @a.orf_length("ATGCCCGGGTAG")
      assert_equal len, 4, "expected 4 but got #{len}"
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

    should "find the mean length" do
      ans = @a.run(2)
      mean = @a.mean_len
      assert_equal mean, 1508.25
    end

  end
end
