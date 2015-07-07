# before the
require 'rbconfig'
require 'yell'
RbConfig::CONFIG['CFLAGS'] = ''

# Transrate is a comprehensive transcriptome assembly
# quality assessment tool.
module Transrate

  # Our own set of errors to allow nice custom error handling
  class TransrateError < StandardError; end
  class TransrateIOError < TransrateError; end
  class TransrateArgError < TransrateError; end

  # Create the universal logger and include it in Object
  # making the logger object available everywhere
  format = Yell::Formatter.new("[%5L] %d : %m", "%Y-%m-%d %H:%M:%S")
  # http://xkcd.com/1179/
  Yell.new(:format => format) do |l|
    l.level = :info
    l.name = Object
    l.adapter STDOUT, level: [:debug, :info, :warn]
    l.adapter STDERR, level: [:error, :fatal]
  end
  Object.send :include, Yell::Loggable

end # Transrate

require 'transrate/transrate'
require 'transrate/transrater'
require 'transrate/version'
require 'transrate/contig'
require 'transrate/assembly'
require 'transrate/snap'
require 'transrate/score_optimiser'
require 'transrate/salmon'
require 'transrate/read_metrics'
require 'transrate/comparative_metrics'
require 'transrate/contig_metrics'
require 'transrate/cmd'
require 'transrate/cmdline'
