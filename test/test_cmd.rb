require 'helper'

class TestCmd < Test::Unit::TestCase

  context "Cmd" do

    should "run commands" do
      cmd = Transrate::Cmd.new 'echo "success"'
      cmd.run
      assert_equal "success", cmd.stdout.chomp, 'run echo command'
    end

  end

end
