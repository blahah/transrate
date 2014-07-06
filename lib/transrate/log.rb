require 'yell'
require 'singleton'

module Transrate

  class Log < Logger

    include Singleton

    def initialize
      logger = Yell.new do |l|
        l.level = :info
        l.adapter STDOUT, level: [:debug, :info, :warn]
        l.adapter STDERR, level: [:error, :fatal]
      end
    end

  end # Log

end # Transrate
