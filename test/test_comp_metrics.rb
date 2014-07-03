require 'helper'

class TestCompMetrics < Test::Unit::TestCase

  context "ComparativeMetrics" do

    setup do
      querypath = File.join(File.dirname(__FILE__),
                            'data',
                            'assembly.2.fa')
      targetpath = File.join(File.dirname(__FILE__),
                            'data',
                            'Os.protein.2.fa')
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
      # Q1   |------------|
      # Q2                    |--------------|
      # T1 |------------------------------------|  # coverage = 200/500
      hash["q1"] << HitHelper.new("q1", "t1", 1, 100, 101, 200, 100, 500)
      hash["q2"] << HitHelper.new("q2", "t1", 1, 100, 301, 400, 100, 500)

      # Q3                    |------------|
      # Q4   |--------------|
      # T2 |------------------------------------| # coverage = 200/500
      hash["q3"] << HitHelper.new("q3", "t2", 1, 100, 301, 400, 100, 500)
      hash["q4"] << HitHelper.new("q4", "t2", 1, 100, 101, 200, 100, 500)

      # Q5                |------------|
      # Q6   |-------------------|
      # T3 |------------------------------------| # coverage = 300/500
      hash["q5"] << HitHelper.new("q5", "t3", 1, 200, 201, 400, 200, 500)
      hash["q6"] << HitHelper.new("q6", "t3", 1, 200, 101, 300, 200, 500)

      # Q7             |------------|
      # Q8      |------------------------|
      # T3 |------------------------------------| # coverage = 300/500
      hash["q7"] << HitHelper.new("q7", "t4", 1, 100, 201, 300, 100, 500)
      hash["q8"] << HitHelper.new("q8", "t4", 1, 300, 101, 400, 300, 500)

      # Q9     |--------|
      # Q10                        |--------|
      # Q11        |--------------------|
      # T5 |------------------------------------| # coverage = 600/1000
      hash["q9"] << HitHelper.new("q9",   "t5", 1, 200, 201, 400, 200, 1000)
      hash["q10"] << HitHelper.new("q10", "t5", 1, 200, 601, 800, 200, 1000)
      hash["q11"] << HitHelper.new("q11", "t5", 1, 400, 301, 700, 400, 1000)

      crb.hash = hash
      ohr = @comp.ortholog_hit_ratio crb
      assert_equal 16.0/30.0, ohr
    end

    should "calculate potential chimera count" do
      crb = CRBHelper.new(false)

      hash = Hash.new
      # (1..3).each do |i|
      #   hash["q#{i}"] = []
      # end
      hash["q1"]=[]

      # T1   |---------|
      # T2                 |---------|
      # Q1 |----------------------------| # chimera = true
      hash["q1"] << HitHelper.new("q1", "t1", 101, 200, 1, 100, 500, 100)
      hash["q1"] << HitHelper.new("q1", "t2", 301, 400, 1, 100, 400, 100)


      # T3   |---------|
      # T3                 |---------|
      # Q2 |----------------------------|
      # chimera = true because the reference has the region 1-100 duplicated
      # hash["q2"] << HitHelper.new("q2", "t3", 101, 200, 1, 100, 500, 100)
      # hash["q2"] << HitHelper.new("q2", "t3", 301, 400, 1, 100, 400, 100)

      # # T3   |---------|
      # # T3                 |---------|
      # # Q2 |----------------------------|
      # # chimera = false because the reference
      # hash["q3"] << HitHelper.new("q3", "t4", 101, 200, 1, 100, 500, 100)
      # hash["q3"] << HitHelper.new("q3", "t4", 301, 400, 1, 100, 400, 100)

      crb.hash = hash
      chi = @comp.chimeras crb
      assert_equal 1, chi
    end

    should "calculate number of contigs with crbblast hit" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir do
          @comp.run
          assert_equal 11, @comp.comp_stats[:n_contigs_with_recip]
          assert_equal 11/13.0, @comp.comp_stats[:p_contigs_with_recip]
        end
      end
    end

    should "calculate number of reference sequences with crbblast hit" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir do
          @comp.run
          assert_equal 10, @comp.comp_stats[:n_refs_with_recip]
          assert_equal 0.5, @comp.comp_stats[:p_refs_with_recip]
        end
      end
    end

    should "calculate reference sequence coverage" do
      # n&p of reference sequences covered to (25, 50, 75, 85, 95%)
      # of their length by CRB-BLAST hit
      Dir.mktmpdir do |tmpdir|
        Dir.chdir do
          @comp.run
          stats = @comp.comp_stats
          assert_equal 10, stats[:cov25]
          assert_equal 10, stats[:cov50]
          assert_equal 7, stats[:cov75]
          assert_equal 6, stats[:cov85]
          assert_equal 3, stats[:cov95]
        end
      end
    end

    should "number of reference sequences coverage" do
      # n&p of reference sequences covered to (25, 50, 75, 85, 95%)
      # of their length by CRB-BLAST hit
      crb = CRBHelper.new(false)

      hash = Hash.new
      (1..5).each do |i|
        hash["q#{i}"] = []
      end
      hash["q1"] << HitHelper.new("q1", "t1", 1, 250, 101, 350, 250, 1000)
      hash["q2"] << HitHelper.new("q2", "t2", 1, 500, 101, 600, 500, 1000)
      hash["q3"] << HitHelper.new("q3", "t3", 1, 750, 101, 850, 750, 1000)
      hash["q4"] << HitHelper.new("q4", "t4", 1, 850, 101, 950, 850, 1000)
      hash["q5"] << HitHelper.new("q5", "t5", 1, 950, 1, 950, 950, 1000)

      crb.hash = hash
      ohr = @comp.ortholog_hit_ratio crb
      stats = @comp.comp_stats
      assert_equal 5, stats[:cov25]
      assert_equal 4, stats[:cov50]
      assert_equal 3, stats[:cov75]
      assert_equal 2, stats[:cov85]
      assert_equal 1, stats[:cov95]
    end

  end
end
