#!/usr/bin/env  ruby

require 'helper'
require 'bio'

class TestInline < Test::Unit::TestCase

  context 'transrate' do

    setup do
      filepath = File.join(File.dirname(__FILE__), 'data', 'assembly.fasta')
      @a = Transrate::Assembly.new(filepath)
    end

    should 'find longest orf in file' do
      orfs = []
      @a.assembly.each_value do |contig|
        orfs << contig.orf_length
      end
      assert_equal 4, orfs.length
      assert_equal [333, 370, 131, 84], orfs
    end

    should 'find longest orf in sequence' do
      seq = Bio::FastaFormat.new ">test\nATGCCCCTAGGGTAG"
      contig = Transrate::Contig.new seq
      assert_equal 4, contig.orf_length
    end

  end
end
