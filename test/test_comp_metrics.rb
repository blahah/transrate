require 'helper'

module CRB_Blast
  class CRB_Blast
    def change_hit(query_name, target_name, qstart, qend, tstart, tend, qlen, tlen)
      hits = @reciprocals[query_name]
      hits.each do |hit|
        if hit.target == target_name
          hit.qstart = qstart
          hit.qend = qend
          hit.tstart = tstart
          hit.tend = tend
          hit.qlen = qlen
          hit.tlen = tlen
        end
      end
    end

    def add_hit(query_name, target_name, qstart, qend, tstart, tend, qlen, tlen)
      @reciprocals[query_name] ||= []
      list = Array.new(14)
      list[0] = query_name
      list[1] = target_name
      list[6] = qstart
      list[7] = qend
      list[8] = tstart
      list[9] = tend
      list[12] = qlen
      list[13] = tlen
      @reciprocals[query_name] << Hit.new(list)
    end

    def remove_hit(query_name)
      @reciprocals.delete(query_name)
    end
  end
end

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

    should "calculate reference coverage" do
      crb = @comp.reciprocal_best_blast
      # change the results so i know what i have
      # qstart, qend, tstart, tend, qlen, tlen
      #
      # Q       |------------|
      # T1 |-------------------------|
      crb.change_hit("scaf_Os03g60760.1", "LOC_Os03g60760.1",
                     1, 300, 101, 200, 300, 200) # 0.5
      @reference["LOC_Os03g60760.1"].seq = "A"*200
      #
      # Q1   |----------|
      # Q2                 |------------|
      # T2  |-------------------------------|
      crb.change_hit("scaf_Os10g39590.1", "LOC_Os10g39590.1",
                     1, 150, 51, 100, 150, 200) # 0.25
      crb.add_hit("scaf_Os10g39590.1", "LOC_Os10g39590.1",
                     1, 150, 151, 200, 150, 200) # 0.25
      @reference["LOC_Os10g39590.1"].seq = "A"*200
      #
      # adding first block [151..300]  scaf_Os09g38670.1
      #    450 / 600.0
      #    LOC_Os09g38670.1  0.75

      #
      #
      # Q1           |-----------|
      # Q2      |----------------------|
      # T3  |-------------------------------|
      crb.change_hit("scaf_Os09g38670.1", "LOC_Os09g38670.1",
                     1, 150, 51, 100, 150, 200) # 0.25
      crb.add_hit("scaf_Os09g38670.1", "LOC_Os09g38670.1",
                     1, 450, 26, 175, 450, 200) # 0.75
      @reference["LOC_Os09g38670.1"].seq = "A"*200

      #
      # Q1      |----------------------|
      # Q2           |-----------|
      # T4  |-------------------------------|
      crb.change_hit("scaf_Os12g21920.1", "LOC_Os12g21920.1", #
                     1, 450, 26, 175, 450, 200) # 0.75
      crb.add_hit("scaf_Os12g21920.1", "LOC_Os12g21920.1",
                     1, 150, 51, 100, 150, 200) # 0.25
      @reference["LOC_Os12g21920.1"].seq = "A"*200


      #
      # Q1    |------|
      # Q2                      |--------|
      # Q3       |-----------------|
      # T5  |-------------------------------|
      crb.change_hit("scaf_Os01g36294.1", "LOC_Os01g36294.1", #
                     1, 300, 51, 100, 300, 400)
      crb.add_hit("scaf_Os01g36294.1", "LOC_Os01g36294.1",
                     1, 300, 200, 250, 300, 400)
      crb.add_hit("scaf_Os01g36294.1", "LOC_Os01g36294.1",
                     1, 300, 75, 225, 300, 400)
      @reference["LOC_Os01g36294.1"].seq = "A"*400

      crb.change_hit("scaf_Os12g22750.1", "LOC_Os12g22750.1",
                     1, 300, 101, 200, 300, 200) # 0.5 # 300/600
      @reference["LOC_Os12g22750.1"].seq = "A"*200

      crb.change_hit("scaf_Os02g55190.1", "LOC_Os02g55190.1",
                     1, 300, 101, 200, 300, 200) # 0.5 # 300/600
      @reference["LOC_Os02g55190.1"].seq = "A"*200

      crb.change_hit("scaf_Os03g56500.1", "LOC_Os03g56500.1",
                     1, 300, 101, 200, 300, 400) # 0.25
      crb.change_hit("scaf_Os03g56500.2", "LOC_Os03g56500.1",
                     1, 300, 201, 300, 300, 400) # 0.25 # 600 / 1200
      @reference["LOC_Os03g56500.1"].seq = "A"*400

      crb.change_hit("scaf_Os03g56724.1", "LOC_Os03g56724.1",
                     1, 300, 101, 200, 300, 200) # 300/600 = 0.5
      @reference["LOC_Os03g56724.1"].seq = "A"*200

      crb.remove_hit("scaf_Os01g11360.1")

      @reference["LOC_Os03g08270.3"].seq = "A"*200
      @reference["LOC_Os10g41970.1"].seq = "A"*200
      @reference["LOC_Os09g26780.1"].seq = "A"*200
      @reference["LOC_Os12g24659.1"].seq = "A"*200
      @reference["LOC_Os01g36410.1"].seq = "A"*200
      @reference["LOC_Os12g22780.1"].seq = "A"*200
      @reference["LOC_Os02g56470.1"].seq = "A"*200
      @reference["LOC_Os03g30530.1"].seq = "A"*200
      @reference["LOC_Os03g49850.1"].seq = "A"*200
      @reference["LOC_Os01g11360.1"].seq = "A"*200
      @reference["LOC_Os01g44140.1"].seq = "A"*200

      assert_equal true, crb.target_is_prot, "target is prot"
      assert_equal false, crb.query_is_prot, "query is prot"
      # total_length of references should be 4400

      cov = @comp.coverage crb
      assert_equal 3600/13200.0, cov, "reference coverage"
    end

    should "calculate potential chimera count" do
      crb = @comp.reciprocal_best_blast
      # # T1   |---------|
      # # T2                 |---------|
      # # Q1 |----------------------------| # chimera = true


      crb.remove_hit("scaf_Os10g39590.1")
      crb.add_hit("scaf_Os10g39590.1", "LOC_Os03g60760.1",
                     1, 150, 51, 100, 400, 60) # 0.25
      crb.add_hit("scaf_Os10g39590.1", "LOC_Os10g39590.1",
                     200, 350, 51, 100, 400, 60) # 0.25

      # # T3   |---------|
      # # T3                 |---------|
      # # Q2 |----------------------------|
      # # chimera = true because the reference has the region 1-100 duplicated
      crb.remove_hit("scaf_Os12g21920.1")
      crb.add_hit("scaf_Os12g21920.1", "LOC_Os12g21920.1",
                     1, 150, 55, 105, 400, 60) # 0.25
      crb.add_hit("scaf_Os12g21920.1", "LOC_Os12g21920.1",
                     200, 350, 51, 100, 400, 60) # 0.25

      @comp.chimeras crb
      assert_equal 2/11.0, @comp.p_chimeras
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
          assert_equal 11, @comp.comp_stats[:n_contigs_with_CRBB]
          assert_equal 11/13.0, @comp.comp_stats[:p_contigs_with_CRBB]
        end
      end
    end

    should "calculate number of reference sequences with crbblast hit" do
      Dir.mktmpdir do |tmpdir|
        Dir.chdir tmpdir do
          @comp.run
          assert_equal 10, @comp.comp_stats[:n_refs_with_CRBB]
          assert_equal 0.5, @comp.comp_stats[:p_refs_with_CRBB]
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

  end
end
