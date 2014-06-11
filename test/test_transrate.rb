#!/usr/bin/env	ruby

require 'helper'

class TestTransrate < Test::Unit::TestCase

  context "transrate" do

    setup do
      filepath = File.join(File.dirname(__FILE__), 'data', 'assembly.fasta')
      @a = Transrate::Assembly.new(filepath)
    end

    should "create assembly object" do
      assert @a
      assert_equal @a.assembly.size, 4
    end

    should "run basic stats" do
      stats = @a.basic_stats
      assert_equal stats["n_seqs"], 4
      assert_equal stats["smallest"], 1409
      assert_equal stats["largest"], 1630
      assert_equal stats["mean_len"], 1508.25
    end

    should "run metrics on assembly" do
      ans = @a.run(2) # using 2 threads
      assert_equal ans, true, "should run but returned #{ans}"
    end

    should "find the mean length" do
      ans = @a.run(2)
      mean = @a.mean_len
      n_bases = @a.n_bases
      assert_equal mean, 1508.25
      assert_equal n_bases, 6033
    end
  end
end
