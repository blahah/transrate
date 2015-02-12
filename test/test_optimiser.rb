require 'helper'
require 'tmpdir'

class TestOptimiser < Test::Unit::TestCase

  context "Optimiser" do

    setup do
      @reference = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.fa')
      left = File.join(File.dirname(__FILE__), "data", "sorghum_100.1.fastq")
      right = File.join(File.dirname(__FILE__), "data", "sorghum_100.2.fastq")
      @assembly = Transrate::Assembly.new @reference
      @readmetrics = Transrate::ReadMetrics.new @assembly

      @readmetrics.run(left, right)
    end

    should "get optimal score" do
      optimiser = Transrate::ScoreOptimiser.new(@assembly, @readmetrics)
      assert_equal 0.1471, optimiser.raw_score.round(4)
      optimal, cutoff = optimiser.optimal_score
      assert_equal 0.4252, optimal.round(4)
      assert_equal 0.5638, cutoff.round(4)
    end

  end
end
