require 'helper'
require 'tmpdir'

class TestReadMetrics < Test::Unit::TestCase

  context "ReadMetrics" do

    setup do
      query = File.join(File.dirname(__FILE__), 'data',
                        'sorghum_transcript.fa')
      assembly = Transrate::Assembly.new(query)
      @read_metrics = Transrate::ReadMetrics.new(assembly)
    end

    teardown do
      Dir.chdir("test") do
        Dir["*bt2"].each do |file|
          File.delete(file)
        end
        Dir["*sam"].each do |file|
          File.delete(file)
        end
      end
    end

    should "check setup ran correctly" do
      assert @read_metrics
    end

    should "run read metrics" do
      left = File.join(File.dirname(__FILE__), 'data', 'left.fastq')
      right = File.join(File.dirname(__FILE__), 'data', 'right.fastq')
      tmpdir = Dir.mktmpdir
      puts tmpdir
      Dir.chdir tmpdir do
        @read_metrics.run(left, right)
        assert @read_metrics.has_run
        stats = @read_metrics.read_stats
        assert_equal 7763, stats[:num_pairs]
      end
    end

    should "find mapping read pairs" do
    end

    should "find read pairs that support the assembly" do
    end

    should "find read pairs that support scaffolding" do
    end

    should "count per-base coverage" do
    end

    should "find median coverage" do
    end

    should "identify contigs with at least one uncovered base" do
    end

    should "identify contigs with coverage <1" do
    end

    should "identify contigs with coverage <10" do
    end

  end

end
