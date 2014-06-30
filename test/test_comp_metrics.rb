#!/usr/bin/env	ruby

require 'helper'

class TestCompMetrics < Test::Unit::TestCase

  context "transrate" do

    setup do
      query = "test/assembly.fasta"
      target = "test/Os.protein.fa"
      assembly = Transrate::Assembly.new(query)
      reference = Transrate::Assembly.new(target)
      threads = 8
      @comp = Transrate::ComparativeMetrics.new(assembly, reference, threads)
    end

    teardown do
      File.delete("test/assembly.nsq") if File.exist?("test/assembly.nsq")
      File.delete("test/assembly.nin") if File.exist?("test/assembly.nin")
      File.delete("test/assembly.nhr") if File.exist?("test/assembly.nhr")
      File.delete("test/Os.protein.psq") if File.exist?("test/Os.protein.psq")
      File.delete("test/Os.protein.pin") if File.exist?("test/Os.protein.pin")
      File.delete("test/Os.protein.phr") if File.exist?("test/Os.protein.phr")
      if File.exist?("test/assembly_into_Os.protein.1.blast")
        File.delete("test/assembly_into_Os.protein.1.blast")
      end
      if File.exist?("test/Os.protein_into_assembly.2.blast")
        File.delete("test/Os.protein_into_assembly.2.blast")
      end
    end

    should "run metrics on assembly" do
      @comp.run
      assert @comp.has_run
    end
  end
end
