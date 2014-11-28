require 'helper'
require 'tmpdir'

class TestExpress < Test::Unit::TestCase

  context "Express" do

    should "load an expression file" do
      file = File.join(File.dirname(__FILE__), 'data',
                        'express_results.xprs')
      e = Transrate::Express.new
      results = e.load_expression file
      assert_equal 4, results.size, "should be four results loaded"
      assert_equal 54, results['C291600'][:eff_len], "eff length is wrong"
      assert_equal 48.005105, results['C291600'][:eff_count],
                   "eff count is wrong"
      assert_equal 5.417487e+00, results['C291600'][:tpm], "tpm is wrong"
    end

  end

end
