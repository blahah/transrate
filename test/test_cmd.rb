require 'helper'

class TestCmd < Test::Unit::TestCase

  context "Cmd" do

    should "run commands" do
      cmd = Transrate::Cmd.new 'echo "success"'
      cmd.run
      assert_equal "success", cmd.stdout.chomp, 'run echo command'
    end

    should "get the string of a command" do
      cmd = Transrate::Cmd.new "echo success"
      assert_equal "echo success", cmd.to_s, 'cmd to string'
    end

  end

end
