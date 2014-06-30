#!/usr/bin/env	ruby

require 'helper'

class TestBowtie < Test::Unit::TestCase

  context "bowtie" do

    setup do
      @reference = "test/sorghum_transcript.fa"
      @left = "test/left.fastq"
      @right = "test/right.fastq"
      @mapper = Transrate::Bowtie2.new

    end

    teardown do
      Dir.chdir("test") do |dir|
        Dir["*bt2"].each do |file|
          File.delete(file)
        end
        Dir["*sam"].each do |file|
          File.delete(file)
        end
      end
    end

    should "build index" do
      @mapper.build_index @reference
      assert File.exist?("test/sorghum_transcript.1.bt2")
    end

    should "build index and map reads" do
      @mapper.build_index @reference
      left = File.basename(@left)
      right = File.basename(@right)
      index = File.basename(@mapper.index_name)
      @mapper.map_reads(@reference, @left, @right)
      sam = @mapper.sam
      assert File.exist?("#{sam}")
    end

    should "raise error when no index built" do
      assert_raise Transrate::Bowtie2Error do
        @mapper.map_reads(@reference, @left, @right)
      end
    end
  end
end