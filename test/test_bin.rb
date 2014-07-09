require 'helper'

class TestTransrateBin < Test::Unit::TestCase

  context "Transrate" do

    setup do

    end

    teardown do
      files = ["150uncovered.l.fq.150uncovered.r.fq.assembly.2.bai",
      "150uncovered.l.fq.150uncovered.r.fq.assembly.2.bam",
      "150uncovered.l.fq.150uncovered.r.fq.assembly.2.sam",
      "150uncovered.l.fq.150uncovered.r.fq.assembly.2.sorted.bam",
      "assembly.2.1.bt2",
      "assembly.2.2.bt2",
      "assembly.2.3.bt2",
      "assembly.2.4.bt2",
      "assembly.2.fa.coverage",
      "assembly.2_into_Os.protein.2.1.blast",
      "assembly.2.nhr",
      "assembly.2.nin",
      "assembly.2.nsq",
      "assembly.2.rev.1.bt2",
      "assembly.2.rev.2.bt2",
      "Os.protein.2_into_assembly.2.2.blast",
      "Os.protein.2.phr",
      "Os.protein.2.pin",
      "Os.protein.2.psq",
      "supported_bridges.csv",
      "transrate.csv"]
      files.each do |file|
        File.delete(file) if File.exist?(file)
      end
    end

    should "run help" do
      c=Transrate::Cmd.new("bundle exec bin/transrate --help")
      c.run
      assert_equal 1755, c.stdout.length, "stdout"
      assert_equal true, c.status.success?, "exit status"
    end

    should "fail on non existent assembly files" do
      c=Transrate::Cmd.new("bundle exec bin/transrate --assembly foo.fasta")
      c.run
      assert_equal 163, c.stderr.length, "stderr"
      assert_equal false, c.status.success?, "exit success"
    end

    should "fail on non existent reference files" do
      c=Transrate::Cmd.new("bundle exec bin/transrate --reference foo.fasta")
      c.run
      assert_equal 104, c.stderr.length, "error"
      assert_equal false, c.status.success?, "exit status"
    end

    should "run on test data" do
      assembly = File.join(File.dirname(__FILE__), 'data', 'assembly.2.fa')
      reference = File.join(File.dirname(__FILE__), 'data', 'Os.protein.2.fa')
      left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      cmd = "bundle exec bin/transrate --assembly #{assembly}"
      cmd << " --reference #{reference}"
      cmd << " --left #{left}"
      cmd << " --right #{right}"
      c = Transrate::Cmd.new("#{cmd}")
      c.run
      assert_equal true, c.status.success?, "exit status"
      assert_equal 1249, File.size("transrate.csv")
    end

    should "fail when one of multiple assemblies is missing" do
      assembly = File.join(File.dirname(__FILE__), 'data', 'assembly.2.fa')
      assembly2 = File.join(File.dirname(__FILE__), 'data', 'foo.fa')
      reference = File.join(File.dirname(__FILE__), 'data', 'Os.protein.2.fa')
      left = File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right = File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      cmd = "bundle exec bin/transrate "
      cmd << " --assembly #{assembly},#{assembly2}"
      cmd << " --reference #{reference}"
      cmd << " --left #{left}"
      cmd << " --right #{right}"
      c = Transrate::Cmd.new("#{cmd}")
      c.run
      assert_equal false, c.status.success?, "exit status"
    end

  end
end
