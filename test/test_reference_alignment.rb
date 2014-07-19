require 'helper'

class TestRefAlignment < Test::Unit::TestCase

  context "ReferenceAlignment" do

    setup do
      querypath = File.join(File.dirname(__FILE__), 'data', 'test_ref_align_assembly.fa')
      targetpath = File.join(File.dirname(__FILE__), 'data', 'test_ref_align_genome.fa')
      assembly = Transrate::Assembly.new(querypath)
      reference = Transrate::Assembly.new(targetpath)
      threads = 4
      @refaln = Transrate::ReferenceAlignment.new(assembly, reference, threads, 1e-5, 97.0, 0)
      @reftest = Transrate::ReferenceAlignment.new(assembly, reference, threads, 1e-5, 99.0, 0)
    end
    should "run metrics on assembly" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @refaln.run
          assert @refaln.has_run
        end
      end
    end
    should "report the number of aligned transcripts" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @refaln.run
          assert_equal 2, @refaln.number_aligned
        end
      end
    end

    should "report the number of unaligned transcripts" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @refaln.run
          assert_equal 1, @refaln.number_unaligned
        end
      end
    end

    should "report the average percentage of aligned bases" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @refaln.run
          assert 96.0 <= @refaln.genome_stats[:avg_perc_aligned]
        end
      end
    end

    should "fail with strict requirement" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @reftest.run
          assert_equal 2, @reftest.number_unaligned
        end
      end
    end
  end
end
