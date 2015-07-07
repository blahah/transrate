require 'helper'
require 'tmpdir'

class TestSnap < MiniTest::Test

  context "Snap" do

    setup do
      @reference = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.fa')
      @left = File.join(File.dirname(__FILE__), "data", "sorghum_100.1.fastq")
      @right = File.join(File.dirname(__FILE__), "data", "sorghum_100.2.fastq")
    end

    should "map reads" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir(tmpdir) do
          snap = Transrate::Snap.new
          snap.build_index(@reference, 4)
          snap.map_reads(@reference, @left, @right)
          assert File.exist?(snap.bam), "bam file not found"
          assert_equal 25006, snap.read_count, "read count"
        end
      end
    end

  end
end
