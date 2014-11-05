# before the
require 'rbconfig'
require 'yell'
RbConfig::CONFIG['CFLAGS'] = ''

require 'transrate/score_optimiser'
require 'transrate/transrater'
require 'transrate/version'
require 'transrate/contig'
require 'transrate/assembly'
require 'transrate/snap'
require 'transrate/express'
require 'transrate/read_metrics'
require 'transrate/comparative_metrics'
require 'transrate/contig_metrics'
require 'transrate/samtools'
require 'transrate/cmd'
require 'transrate/sam_checker'
require 'transrate/transrate.so'

# Transrate is a comprehensive transcriptome assembly
# quality assessment tool.
module Transrate

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
