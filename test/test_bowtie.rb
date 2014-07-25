#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestBowtie < Test::Unit::TestCase

  context "bowtie" do

    setup do
      @reference = File.join(File.dirname(__FILE__), 'data',
       'sorghum_transcript.fa')
      @left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      @right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      @unpaired = File.join(File.dirname(__FILE__), 'data', '150uncovered.u.fq')
	  @library = "fr"
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
          unpaired = File.basename(@unpaired)
          index = File.basename(@mapper.index_name)
          @mapper.map_reads(@reference, @left, @right, @unpaired, @library)
          sam = @mapper.sam
          assert File.exist?("#{sam}"), "sam file doesn't exist"
          cmd = "grep -v \"^@\" #{sam} | wc -l "
          line_in_sam_file = `#{cmd}`.chomp.to_i
          assert_equal 471, line_in_sam_file, "lines in sam file"
        end
      end
    end

    should "raise error when no index built" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          assert_raise Transrate::Bowtie2Error do
            @mapper.map_reads(@reference, @left, @right, @unpaired, @library)
          end
        end
      end
    end

    should "raise error when bowtie fails" do
      not_reads = File.join(File.dirname(__FILE__), 'data', 'not_a_file.fq')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          assert_raise Transrate::Bowtie2Error do
            @mapper.build_index @reference
            @mapper.map_reads(@reference, @left, not_reads, nil, @library)
          end
        end
      end
    end
  end
end
