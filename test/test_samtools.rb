require 'helper'

class TestSamtools < Test::Unit::TestCase

  context "samtools" do

    should "know the path to samtools binary" do
      msg = /Program: samtools/
      path = Transrate::Samtools.path
      res = `#{path} 2>&1`.split("\n").join
      assert msg =~ res
    end

    should "run commands" do
      sam = File.join(File.dirname(__FILE__), 'data', 'tiny.sam')
      Transrate::Samtools.run "view -bS #{sam} > tiny.bam"
      assert_equal 460, File.size('tiny.bam'), 'bam file should be created'
      File.delete 'tiny.bam'
    end

  end
end
