require 'helper'

class TestCompMetrics < Test::Unit::TestCase

  context "transrate" do

    setup do
      querypath = File.join(File.dirname(__FILE__),
                            'data',
                            'assembly.fa')
      targetpath = File.join(File.dirname(__FILE__), 'data', 'Os.protein.fa')
      assembly = Transrate::Assembly.new(querypath)
      reference = Transrate::Assembly.new(targetpath)
      threads = 8
      @comp = Transrate::ComparativeMetrics.new(assembly, reference, threads)
    end

    teardown do
      File.delete("test/assembly.nsq") if File.exist?("test/assembly.nsq")
      File.delete("test/assembly.nin") if File.exist?("test/assembly.nin")
      File.delete("test/assembly.nhr") if File.exist?("test/assembly.nhr")
      File.delete("test/Os.protein.psq") if File.exist?("test/Os.protein.psq")
      File.delete("test/Os.protein.pin") if File.exist?("test/Os.protein.psq")
      File.delete("test/Os.protein.phr") if File.exist?("test/Os.protein.psq")
      if File.exist?("assembly_into_Os.protein.1.blast")
        File.delete("assembly_into_Os.protein.1.blast")
      end
      if File.exist?("Os.protein_into_assembly.2.blast")
        File.delete("Os.protein_into_assembly.2.blast")
      end
    end

    should "run metrics on assembly" do
      @comp.run
      assert @comp.has_run
    end
  end
end
