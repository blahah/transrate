require 'helper'
require 'csv'

class TestCmdline < MiniTest::Test

  context "Transrate" do

    teardown do
      `rm -rf transrate_results`
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
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          assembly, reference, left, right = sorghum_data
          assert_raises Transrate::TransrateIOError do
            Transrate::Cmdline.new("--assembly #{assembly} --reference foo.fasta".split)
          end
        end
      end
    end

    should "run on test data" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          captured_stdout = capture_stdout do
            assembly, reference, left, right = sorghum_data
            cmd = "--assembly #{assembly}"
            cmd << " --reference #{reference}"
            cmd << " --left #{left}"
            cmd << " --right #{right}"
            c = Transrate::Cmdline.new(cmd.split)
            c.run
          end
          assert File.exist?("transrate_results/assemblies.csv"), "csv file doesn't exist"
          assert File.exist?("transrate_results/sorghum_100/contigs.csv"),
                 "contig csv file doesn't exist"
          hash = {}
          CSV.foreach("transrate_results/assemblies.csv", :headers => true,
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
      end
    end

    should "run on test data with comma separated list of fastq files" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          captured_stdout = capture_stdout do
            assembly, reference, left, right = sorghum_data
            left = [left]
            left << File.join(File.dirname(__FILE__), 'data', '150uncovered.l.fq')
            right = [right]
            right << File.join(File.dirname(__FILE__), 'data', '150uncovered.r.fq')
            cmd = "--assembly #{assembly}"
            cmd << " --left #{left.join(",")}"
            cmd << " --right #{right.join(",")}"
            c = Transrate::Cmdline.new(cmd.split)
            c.run
          end
          assert File.exist?("transrate_results/assemblies.csv")
          assert File.exist?("transrate_results/sorghum_100/contigs.csv"),
                 "contig csv file doesn't exist"
          hash = {}
          CSV.foreach("transrate_results/assemblies.csv", :headers => true,
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
      end
    end

    should "fail when one of multiple assemblies is missing" do
      assembly, reference, left, right = sorghum_data
      assembly2 = "does_not_exist.fa"
      cmd = " --assembly #{assembly},#{assembly2}"
      cmd << " --reference #{reference}"
      cmd << " --left #{left}"
      cmd << " --right #{right}"
      assert_raises Transrate::TransrateIOError do
        Transrate::Cmdline.new(cmd.split)
      end
    end

    should "fail when assemblies output already exists" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          Dir.mkdir("test")
          File.open("test/assemblies.csv", "w") { |out| out.write "test\n" }
          assembly, reference, left, right = sorghum_data
          cmd = " --assembly #{assembly}"
          cmd << " --reference #{reference}"
          cmd << " --left #{left}"
          cmd << " --right #{right}"
          cmd << " --output test"
          assert_raises Transrate::TransrateArgError do
            Transrate::Cmdline.new(cmd.split)
          end
        end
      end
    end

  end
end
