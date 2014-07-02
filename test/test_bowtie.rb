#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestBowtie < Test::Unit::TestCase

  context "bowtie" do

    setup do
      @reference = File.join(File.dirname(__FILE__), 'data',
       'sorghum_transcript.fa')
      @left = File.join(File.dirname(__FILE__), 'data', 'left.fastq')
      @right = File.join(File.dirname(__FILE__), 'data', 'right.fastq')
      @mapper = Transrate::Bowtie2.new
    end

    should "build index" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @mapper.build_index @reference
          assert File.exist?("sorghum_transcript.1.bt2")
        end
      end
    end

    should "build index and map reads" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @mapper.build_index @reference
          left = File.basename(@left)
          right = File.basename(@right)
          index = File.basename(@mapper.index_name)
          @mapper.map_reads(@reference, @left, @right)
          sam = @mapper.sam
          assert File.exist?("#{sam}")
          assert_equal 4940290, File.size(sam)
        end
      end
    end

    should "raise error when no index built" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          assert_raise Transrate::Bowtie2Error do
            @mapper.map_reads(@reference, @left, @right)
          end
        end
      end
    end
  end
end
