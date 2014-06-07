#!/usr/bin/env	ruby

require 'helper'

class TestReciprocalAnnotation < Test::Unit::TestCase

  context "reciprocal annotation" do

    setup do
      querypath = File.join(File.dirname(__FILE__),
                            'data',
                            'chimeric_contig.fa')
      targetpath = File.join(File.dirname(__FILE__), 'data', 'at_p.faa')
      assembly = Transrate::Assembly.new(querypath)
      reference = Transrate::Assembly.new(targetpath)
      @recip = Transrate::ReciprocalAnnotation.new(assembly, reference)
    end

    should "run reciprocal annotation metrics on assembly" do
      assert @recip.run
    end

  end
end
