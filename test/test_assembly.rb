require 'helper'
require 'bio'

class TestAssembly < Test::Unit::TestCase

  context "Assembly" do

    setup do
      a = File.join(File.dirname(__FILE__), 'data', 'sorghum_transcript.fa')
      @assembly = Transrate::Assembly.new(a)
    end

    should "run" do
      @assembly.run 1
      assert_equal true, @assembly.has_run
    end

    should "classify contigs" do
      @assembly.run 1
      read_metrics = Transrate::ReadMetrics.new(@assembly)
      left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          read_metrics.run(left, right)
          @assembly.classify_contigs
          assert_equal 840, File.stat("good.sorghum_transcript.fa").size
          assert_equal 0, File.stat("fragmented.sorghum_transcript.fa").size
          assert_equal 749, File.stat("chimeric.sorghum_transcript.fa").size
          assert_equal 0, File.stat("bad.sorghum_transcript.fa").size
        end
      end
    end

  end
end
