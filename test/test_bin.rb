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
      "test_ref_align_assembly.nhr",
      "test_ref_align_assembly.nin",
      "test_ref_align_assembly.nsq",
      "test_ref_align_assembly_into_test_ref_align_genome.1.blast",
      "test_ref_align_genome.nhr",
      "test_ref_align_genome.nin",
      "test_ref_align_genome.nsq",
      "test_ref_align_genome_into_test_ref_align_assembly.2.blast",
      "Os.protein.2_into_assembly.2.2.blast",
      "Os.protein.2.phr",
      "Os.protein.2.pin",
      "Os.protein.2.psq",
      "supported_bridges.csv",
      "singletons.fa",
      "newfile.sam",
      "unaligned_transcripts.fa",
      "transrate.csv"]
      files.each do |file|
        File.delete(file) if File.exist?(file)
      end
    end

    should "run help" do
      c=Transrate::Cmd.new("bundle exec bin/transrate --help")
      c.run
      assert_equal 2315, c.stdout.length, "stdout"
      assert_equal true, c.status.success?, "exit status"
    end

    should "fail on non existent assembly files" do
      c=Transrate::Cmd.new("bundle exec bin/transrate --assembly foo.fasta")
      c.run
      assert_equal 166, c.stderr.length, "stderr"
      assert_equal false, c.status.success?, "exit success"
    end

    should "fail on non existent reference files" do
      c=Transrate::Cmd.new("bundle exec bin/transrate --reference foo.fasta")
      c.run
      assert_equal 105, c.stderr.length, "error"
      assert_equal false, c.status.success?, "exit status"
    end

    should "fail on non existent genome files" do
      assembly = File.join(File.dirname(__FILE__), 'data', 'assembly.2.fa')
      c=Transrate::Cmd.new("bundle exec bin/transrate --assembly #{assembly} --genome foo.fasta")
      c.run
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
      cmd << " --singletons newfile.sam"
      c = Transrate::Cmd.new("#{cmd}")
      c.run
      assert_equal true, c.status.success?, "exit status"
      assert File.exist?("transrate.csv")
      hash = {}
      CSV.foreach("transrate.csv", :headers => true,
                                   :header_converters => :symbol,
                                   :converters => :all) do |row|
        row.headers
        row.fields
        row.headers.zip(row.fields).each do |header, field|
          hash[header]=field
        end
      end
      assert_equal 10331, hash[:n_bases]
      assert_equal 1566, hash[:n50]
      assert_equal 10, hash[:n_refs_with_recip]
      assert_equal 11, `wc -l newfile.sam`.split[0].to_i - 4
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

    should "run with the genome argument provided" do
      assembly = File.join(File.dirname(__FILE__), 'data', 'test_ref_align_assembly.fa')
      genome = File.join(File.dirname(__FILE__), 'data', 'test_ref_align_genome.fa')
      cmd = "bundle exec bin/transrate "
      cmd << " --assembly #{assembly}"
      cmd << " --genome #{genome}"
      c = Transrate::Cmd.new("#{cmd}")
      c.run
      assert_equal true, c.status.success?, "exit status"
    end

    should "report on the assembly's alignment to the genome" do
      assembly = File.join(File.dirname(__FILE__), 'data', 'test_ref_align_assembly.fa')
      genome = File.join(File.dirname(__FILE__), 'data', 'test_ref_align_genome.fa')
      cmd = "bundle exec bin/transrate "
      cmd << " --assembly #{assembly}"
      cmd << " --genome #{genome}"
      c = Transrate::Cmd.new("#{cmd}")
      c.run
      assert File.exist?("transrate.csv")
      hash={}
      CSV.foreach("transrate.csv", :headers => true,
                                   :header_converters => :symbol,
                                   :converters => :all) do |row|
        row.headers
        row.fields
        row.headers.zip(row.fields).each do |header, field|
          hash[header]=field
        end
      end
      assert_equal 3, hash[:aligned]
    end

  end
end
