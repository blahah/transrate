require 'helper'

class TestCompMetrics < Test::Unit::TestCase

  context "ComparativeMetrics" do

    setup do
      querypath = File.join(File.dirname(__FILE__),
                            'data',
                            'assembly.fasta')
      targetpath = File.join(File.dirname(__FILE__),
                            'data',
                            'Os.protein.fa')
      assembly = Transrate::Assembly.new(querypath)
      reference = Transrate::Assembly.new(targetpath)
      threads = 8
      @comp = Transrate::ComparativeMetrics.new(assembly, reference, threads)
    end


    should "run metrics on assembly" do
      @comp.run
      assert @comp.has_run
    end
  end
end
