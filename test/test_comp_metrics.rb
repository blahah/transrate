#!/usr/bin/env	ruby

require 'helper'

class TestCompMetrics < Test::Unit::TestCase

  context "transrate" do

    setup do
      # @a = Transrate::Assembly.new("test/assembly.fasta")
      query = "test/Alyrata.fa"
      target = "test/Athaliana.fa"
      assembly = Transrate::Assembly.new(query)
      reference = Transrate::Assembly.new(target)
      @comp = Transrate::ComparativeMetrics.new(assembly, reference)
    end

    should "run metrics on assembly" do
      @comp.run
      assert @comp.has_run
    end

  end
end
