#!/usr/bin/env	ruby

require 'helper'

class TestReciprocalAnnotation < Test::Unit::TestCase

  context "reciprocal annotation" do

    setup do
      # @a = Transrate::Assembly.new("test/assembly.fasta")
      query = "test/chimeric_contig.fa"
      target = "test/at_p.faa"
      assembly = Transrate::Assembly.new(query)
      reference = Transrate::Assembly.new(target)
      @recip = Transrate::ReciprocalAnnotation.new(assembly, reference)
    end

    should "run reciprocal annotation metrics on assembly" do
      assert @recip.run
    end

  end
end
