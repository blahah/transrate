require 'logger'

module Transrate

  class Log < Logger

    def dump_process?(status, output, msg)
      unless status.exitstatus == 0
        fatal(msg)
        fatal(output)
      end
    end
    
  end # Log

end # Transrate
