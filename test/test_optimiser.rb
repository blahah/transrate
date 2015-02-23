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
      assert_in_delta 0.1507, optimiser.raw_score.round(4), 0.005
      optimal, cutoff = optimiser.optimal_score
      assert_in_delta 0.4529, optimal.round(4), 0.05
      assert_in_delta 0.5638, cutoff.round(4), 0.04
    end

  end
end
