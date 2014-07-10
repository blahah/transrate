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
      @assembly = Transrate::Assembly.new(querypath)
      @q_ids = @assembly.assembly.keys
      @reference = Transrate::Assembly.new(targetpath)
      @t_ids = @reference.assembly.keys
      threads = 8
      @comp = Transrate::ComparativeMetrics.new(@assembly, @reference, threads)
    end


    should "run metrics on assembly" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @comp.run
          assert @comp.has_run
        end
      end
    end

    should "calculate ortholog hit ratio" do
      crb = CRBHelper.new(true)
      hash = Hash.new([])

      # Q1   |------------|
      # Q2                    |--------------|
      # T1 |------------------------------------|  # coverage = 200/500
      hash[@q_ids[0]] <<
        HitHelper.new(@q_ids[0], @t_ids[0], 1, 100, 101, 133, 100, 500)
      hash[@q_ids[1]] <<
        HitHelper.new(@q_ids[1], @t_ids[0], 1, 100, 301, 333, 100, 500)

      # Q3                    |------------|
      # Q4   |--------------|
      # T2 |------------------------------------| # coverage = 200/500
      hash[@q_ids[2]] <<
        HitHelper.new(@q_ids[2], @t_ids[1], 1, 100, 301, 433, 100, 500)
      hash[@q_ids[3]] <<
        HitHelper.new(@q_ids[4], @t_ids[1], 1, 100, 101, 233, 100, 500)

      # Q5                |------------|
      # Q6   |-------------------|
      # T3 |------------------------------------| # coverage = 300/500
      hash[@q_ids[4]] <<
        HitHelper.new(@q_ids[4], @t_ids[2], 1, 200, 201, 233, 200, 500)
      hash[@q_ids[5]] <<
        HitHelper.new(@q_ids[5], @t_ids[2], 1, 200, 101, 133, 200, 500)

      # Q7             |------------|
      # Q8      |------------------------|
      # T3 |------------------------------------| # coverage = 300/500
      hash[@q_ids[6]] <<
        HitHelper.new(@q_ids[6], @t_ids[3], 1, 100, 201, 233, 100, 500)
      hash[@q_ids[7]] <<
        HitHelper.new(@q_ids[7], @t_ids[3], 1, 300, 101, 133, 300, 500)

      # Q9     |--------|
      # Q10                        |--------|
      # Q11        |--------------------|
      # T5 |------------------------------------| # coverage = 600/1000
      hash[@q_ids[8]] <<
        HitHelper.new(@q_ids[8],   @t_ids[4], 1, 200, 201, 400, 200, 1000)
      hash[@q_ids[9]] <<
        HitHelper.new(@q_ids[9], @t_ids[4], 1, 200, 601, 800, 200, 1000)
      hash[@q_ids[10]] <<
        HitHelper.new(@q_ids[10], @t_ids[4], 1, 400, 301, 700, 400, 1000)

      crb.hash = hash

      ohr = @comp.ortholog_hit_ratio crb
      assert_equal 16.0/30.0, ohr
    end

    should "calculate potential chimera count" do
      crb = CRBHelper.new(true)
      hash = Hash.new([])

      # T1   |---------|
      # T2                 |---------|
      # Q1 |----------------------------| # chimera = true
      hash[@q_ids[0]] <<
        HitHelper.new(@q_ids[0], @t_ids[0], 101, 200, 1, 100, 500, 100)
      hash[@q_ids[0]] <<
        HitHelper.new(@q_ids[0], @t_ids[1], 301, 400, 1, 100, 400, 100)


      # T3   |---------|
      # T3                 |---------|
      # Q2 |----------------------------|
      # chimera = true because the reference has the region 1-100 duplicated
      hash[@q_ids[1]] <<
        HitHelper.new(@q_ids[1], @t_ids[2], 101, 200, 1, 100, 500, 100)
      hash[@q_ids[1]] <<
        HitHelper.new(@q_ids[1], @t_ids[2], 301, 400, 1, 100, 400, 100)

      # # T4   |---------|
      # # T4                 |---------|
      # # Q3 |----------------------------|
      # # chimera = false because the reference
      hash[@q_ids[2]] <<
        HitHelper.new(@q_ids[2],@t_ids[3], 101, 200, 1, 100, 500, 200)
      hash[@q_ids[2]] <<
        HitHelper.new(@q_ids[2], @t_ids[3], 301, 400, 101, 200, 400, 200)

      crb.hash = hash

      @comp.chimeras crb
      assert_equal 0.667, @comp.p_chimeras.round(3)
    end

    should "calculate overlap amount" do
      assert_equal 0.5, @comp.overlap_amount(201,500,101,400), "1"
      assert_equal 0.5, @comp.overlap_amount(101,400,201,500), "2"
      assert_equal 0.5, @comp.overlap_amount(201,400,101,500), "3"
      assert_equal 0.5, @comp.overlap_amount(101,500,201,400), "4"
    end

    should "calculate number of contigs with crbblast hit" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @comp.run
          assert_equal 11, @comp.comp_stats[:n_contigs_with_recip]
          assert_equal 11/13.0, @comp.comp_stats[:p_contigs_with_recip]
        end
      end
    end

    should "calculate number of reference sequences with crbblast hit" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
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
        Dir.chdir tmpdir do
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
      crb = CRBHelper.new(true)
  hash = Hash.new([])

      hash[@q_ids[0]] <<
        HitHelper.new(@q_ids[0], @t_ids[0], 1, 250, 101, 350, 250, 1000)
      hash[@q_ids[1]] <<
        HitHelper.new(@q_ids[1], @t_ids[1], 1, 500, 101, 600, 500, 1000)
      hash[@q_ids[2]] <<
        HitHelper.new(@q_ids[2], @t_ids[2], 1, 750, 101, 850, 750, 1000)
      hash[@q_ids[3]] <<
        HitHelper.new(@q_ids[3], @t_ids[3], 1, 850, 101, 950, 850, 1000)
      hash[@q_ids[4]] <<
        HitHelper.new(@q_ids[4], @t_ids[4], 1, 950, 1, 950, 950, 1000)

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
