#!/usr/bin/env	ruby

require 'helper'

class TestTransrate < Test::Unit::TestCase

  context "transrate" do

    setup do
      filepath = File.join(File.dirname(__FILE__), 'data', 'assembly.fasta')
      @a = Transrate::Assembly.new(filepath)
    end

    should "run metrics on assembly" do
      ans = @a.run(2) # using 2 threads
      assert_equal ans, true, "should run but returned #{ans}"
    end

    should "find the mean length" do
      ans = @a.run(2)
      mean = @a.mean_len
      assert_equal mean, 1508.25
    end
  end
end
