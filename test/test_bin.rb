require 'helper'
require 'csv'

class TestTransrateBin < Test::Unit::TestCase

  context "Transrate" do

    setup do

    end

    teardown do
      files = ["sorghum_100.1.fastq.sorghum_100.2.fastq.sorghum_100.bam.bai",
      "sorghum_100.1.fastq.sorghum_100.2.fastq.sorghum_100.bam",
      "sorghum_100_into_Os.protein.2.1.blast",
      "assembly.2_into_Os.protein.2.1.blast",
      "sorghum_100.nhr", "sorghum_100.nin", "sorghum_100.nsq",
      "Os.protein.2_into_sorghum_100.2.blast",
      "Os.protein.2_into_assembly.2.2.blast",
      "assembly.2.nhr",  "assembly.2.nin",  "assembly.2.nsq",
      "Os.protein.2.phr",  "Os.protein.2.pin",  "Os.protein.2.psq",
      "transrate_assemblies.csv",
      "sorghum_100.fa_quant.sf",
      "sorghum_100.fa_bam_info.csv",
      "transrate_sorghum_100.fa_contigs.csv",
      "sorghum_100.1.fastq-sorghum_100.2.fastq-read_count.txt",
      "sorghum_100.1.fastq.sorghum_100.2.fastq.sorghum_100.assigned.bam",
      "bad.sorghum_100.fa", "chimeric.sorghum_100.fa",
      "fragmented.sorghum_100.fa", "good.sorghum_100.fa"]
      files.each do |file|
        File.delete(file) if File.exist?(file)
      end
      `rm -rf sorghum_100`
      `rm -rf logs`
      `rm -rf libParams`
    end

    should "run help" do
      captured_stdout = capture_stdout do
        Transrate::Cmdline.new("--help".split)
      end
      at_exit do
        assert last_exit_successful?, "exit success"
        assert captured_stdout =~ /Analyse a de-novo transcriptome assembly/,
               "error message"
      end
    end

    should "fail nicely on non existent assembly files" do
      assert_raises Transrate::TransrateIOError do
        Transrate::Cmdline.new("--assembly foo.fasta".split)
      end
    end

    should "fail nicely when assembly is not provided" do
      assert_raises Transrate::TransrateArgError do
        Transrate::Cmdline.new("--left left.fq".split)
      end
    end

    should "fail on non existent reference files" do
      assembly, reference, left, right = sorghum_data
      assert_raises Transrate::TransrateIOError do
        Transrate::Cmdline.new("--assembly #{assembly} --reference foo.fasta".split)
      end
    end

    should "run on test data" do
      assembly, reference, left, right = sorghum_data
      cmd = "--assembly #{assembly}"
      cmd << " --reference #{reference}"
      cmd << " --left #{left}"
      cmd << " --right #{right}"
      c = Transrate::Cmdline.new(cmd.split)
      c.run
      at_exit do

      end
      assert_equal true, last_exit_successful?, "exit status"
      assert File.exist?("transrate_assemblies.csv"), "csv file doesn't exist"
      assert File.exist?("transrate_sorghum_100.fa  _contigs.csv"),
             "contig csv file doesn't exist"
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
      assert_in_delta 137748, hash[:n_bases], 1000, "number of bases"
      assert_equal 1692, hash[:n50], "n50"
      assert_equal 2, hash[:n_refs_with_crbb], "number of crb hits"
      assert_equal 2, hash[:n_contigs_with_crbb], "number of contigs with hits"
    end

    should "run on test data with comma separated list of fastq files" do
      assembly, reference, left, right = sorghum_data
      left = [left]
      left << File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
      right = [right]
      right << File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
      cmd = "--assembly #{assembly}"
      cmd << " --left #{left.join(",")}"
      cmd << " --right #{right.join(",")}"
      c = Transrate::Cmdline.new("#{cmd}")
      c.run
      if !(c.status.success?)
        puts c.stdout
        puts c.stderr
      end
      assert_equal true, c.status.success?, "exit status"
      assert File.exist?("transrate_assemblies.csv")
      assert File.exist?("transrate_sorghum_100.fa_contigs.csv"),
             "contig csv file doesn't exist"
      hash = {}
      CSV.foreach("transrate_assemblies.csv", :headers => true,
                                   :header_converters => :symbol,
                                   :converters => :all) do |row|
        row.headers.zip(row.fields).each do |header, field|
          hash[header]=field
        end
      end
      assert_in_delta 137748, hash[:n_bases], 1000, "number of bases"
      assert_equal 1692, hash[:n50], "n50"
      assert_equal 25006 + 223, hash[:fragments], "number of reads"
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
