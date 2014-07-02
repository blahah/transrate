require 'helper'

class TestCompMetrics < Test::Unit::TestCase

  context "ComparativeMetrics" do

    setup do
      querypath = File.join(File.dirname(__FILE__),
                            'data',
                            'assembly.fasta')
      targetpath = File.join(File.dirname(__FILE__),
                            'data',
                            'Os.protein.fa')
      assembly = Transrate::Assembly.new(querypath)
      reference = Transrate::Assembly.new(targetpath)
      threads = 8
      @comp = Transrate::ComparativeMetrics.new(assembly, reference, threads)
    end


    should "run metrics on assembly" do
      @comp.run
      assert @comp.has_run
    end

    should "calculate ortholog hit ratio" do
      crb = CRBHelper.new(false)

      hash = Hash.new
      (1..11).each do |i|
        hash["q#{i}"] = []
      end
      hash["q1"] << HitHelper.new("q1", "t1", 101, 300, 1, 400, 400, 400)
      hash["q2"] << HitHelper.new("q2", "t1", 201, 400, 1, 400, 400, 400)
      hash["q3"] << HitHelper.new("q3", "t2", 110, 290, 1, 400, 400, 400)
      hash["q4"] << HitHelper.new("q4", "t2", 101, 300, 1, 400, 400, 400)
      hash["q5"] << HitHelper.new("q5", "t3", 200, 400, 1, 400, 400, 400)
      hash["q6"] << HitHelper.new("q6", "t3", 101, 300, 1, 400, 400, 400)
      hash["q7"] << HitHelper.new("q7", "t4", 101, 200, 1, 400, 400, 400)
      hash["q8"] << HitHelper.new("q8", "t4", 301, 400, 1, 400, 400, 400)
      hash["q9"] << HitHelper.new("q9", "t5", 101, 200, 1, 400, 400, 400)
      hash["q10"] << HitHelper.new("q10", "t5", 301, 400, 1, 400, 400, 400)
      hash["q11"] << HitHelper.new("q11", "t5", 150, 350, 1, 400, 400, 400)
      crb.hash = hash
      ohr = @comp.ortholog_hit_ratio crb
      assert_equal 0.65, ohr
    end

    should "calculate potential chimera count" do
    end

    should "calculate number of contigs with crbblast hit" do
    end

    should "calculate number of reference sequences with crbblast hit" do
    end

    should "calculate references sequences coverage" do
    end

  end
end
