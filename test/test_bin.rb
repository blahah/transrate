require 'helper'
require 'csv'

class TestTransrateBin < Test::Unit::TestCase

  context "Transrate" do

    setup do

    end

    teardown do
      files = ["150uncovered.l.fq.150uncovered.r.fq.assembly.2.bai",
      "150uncovered.l.fq.150uncovered.r.fq.assembly.2.bam",
      "150uncovered.l.fq.150uncovered.r.fq.assembly.2.sam",
      "150uncovered.l.fq.150uncovered.r.fq.assembly.2.sorted.bam",
      "assembly.2.1.bt2", "assembly.2.2.bt2", "assembly.2.3.bt2",
      "assembly.2.4.bt2", "assembly.2.fa.coverage",
      "assembly.2_into_Os.protein.2.1.blast",
      "assembly.2.nhr", "assembly.2.nin", "assembly.2.nsq",
      "assembly.2.rev.1.bt2",  "assembly.2.rev.2.bt2",
      "Os.protein.2_into_assembly.2.2.blast",
      "Os.protein.2.phr",  "Os.protein.2.pin",  "Os.protein.2.psq",
      "supported_bridges.csv",
      "transrate_assemblies.csv",
      "transrate_contigs.csv","assembly.2.fa.bcf",
      "transrate_assembly.2.fa_contigs.csv",
      "test.sam"]
      files.each do |file|
        File.delete(file) if File.exist?(file)
      end
    end

    should "run help" do
      c=Transrate::Cmd.new("bundle exec bin/transrate --help")
      c.run
      assert_equal 1818, c.stdout.length, "stdout"
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
      cmd << " --singletons test.sam"
      c = Transrate::Cmd.new("#{cmd}")
      c.run
      assert_equal true, c.status.success?, "exit status"
      assert File.exist?("transrate_assemblies.csv"), "csv file doesn't exist"
      assert File.exist?("transrate_assembly.2.fa_contigs.csv"),
             "csv file doesn't exist"
      hash = {}
      CSV.foreach("transrate_assemblies.csv", :headers => true,
                                   :header_converters => :symbol,
                                   :converters => :all) do |row|
        row.headers
        row.fields
        row.headers.zip(row.fields).each do |header, field|
          hash[header]=field
        end
      end
      assert_equal 10331, hash[:n_bases], "number of bases"
      assert_equal 1566, hash[:n50], "n50"
      assert_equal 10, hash[:n_refs_with_crbb], "number of crb hits"
      assert_equal 11, `wc -l test.sam`.split[0].to_i - 4
    end

    should "run on test data with comma separated list of fastq files" do
      assembly = File.join(File.dirname(__FILE__), 'data', 'assembly.2.fa')
      left = []
      right = []
      left << File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      left << File.join(File.dirname(__FILE__), 'data', 'bridging_reads.l.fastq')
      right << File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      right << File.join(File.dirname(__FILE__), 'data', 'bridging_reads.r.fastq')
      cmd = "bundle exec bin/transrate --assembly #{assembly}"
      cmd << " --left #{left.join(",")}"
      cmd << " --right #{right.join(",")}"
      c = Transrate::Cmd.new("#{cmd}")
      c.run
      assert_equal true, c.status.success?, "exit status"
      assert File.exist?("transrate_assemblies.csv"), "csv file doesn't exist"
      assert File.exist?("transrate_assembly.2.fa_contigs.csv"),
             "csv file doesn't exist"
      hash = {}
      CSV.foreach("transrate_assemblies.csv", :headers => true,
                                   :header_converters => :symbol,
                                   :converters => :all) do |row|
        row.headers
        row.fields
        row.headers.zip(row.fields).each do |header, field|
          hash[header]=field
        end
      end
      assert_equal 10331, hash[:n_bases], "number of bases"
      assert_equal 1566, hash[:n50], "n50"
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
