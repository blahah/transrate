#!/usr/bin/env	ruby

require 'helper'

class TestReadMetrics < Test::Unit::TestCase

  context "transrate" do

    setup do
      query = "test/sorghum_transcript.fa"
      assembly = Transrate::Assembly.new(query)
      @read_metrics = Transrate::ReadMetrics.new(assembly)
    end

    teardown do
      Dir.chdir("test") do
        Dir["*bt2"].each do |file|
          File.delete(file)
        end
        Dir["*sam"].each do |file|
          #File.delete(file)
        end
      end
    end

    should "check setup ran correctly" do
      assert @read_metrics
    end

    should "run read metrics" do
      left = "test/left.fastq"
      right = "test/right.fastq"
      @read_metrics.run(left, right)
      assert @read_metrics.has_run
      stats = @read_metrics.read_stats
      assert_equal 7763, stats[:num_pairs]
      p stats
    end

  end

end