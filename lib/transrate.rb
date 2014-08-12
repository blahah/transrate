# before the
require 'rbconfig'
require 'yell'
RbConfig::CONFIG['CFLAGS'] = ''

require 'transrate/transrater'
require 'transrate/version'
require 'transrate/contig'
require 'transrate/assembly'
require 'transrate/bowtie2'
require 'transrate/read_metrics'
require 'transrate/comparative_metrics'
require 'transrate/contig_metrics'
require 'transrate/metric'
require 'transrate/dimension_reduce'
require 'transrate/samtools'
require 'transrate/cmd'
require 'transrate/transrate.so'

# Transrate is a comprehensive transcriptome assembly
# quality assessment tool.
module Transrate

  # Create the universal logger and include it in Object
  # making the logger object available everywhere
  Yell.new(:format => "[%5L]: %m") do |l|
    l.level = :info
    l.name = Object
    l.adapter STDOUT, level: [:debug, :info, :warn]
    l.adapter STDERR, level: [:error, :fatal]
  end
  Object.send :include, Yell::Loggable

end # Transrate
