#!/usr/bin/env	ruby

require 'helper'

class TestCompMetrics < Test::Unit::TestCase

  context "transrate" do

    setup do
      query = "test/clg.head.assembly.fa"
      target = "test/Os.protein.fa"
      assembly = Transrate::Assembly.new(query)
      reference = Transrate::Assembly.new(target)
      threads = 8
      @comp = Transrate::ComparativeMetrics.new(assembly, reference, threads)
    end

    should "run metrics on assembly" do
      @comp.run
      assert @comp.has_run
    end
  end
end
