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

# Transrate is a comprehensive transcriptome assembly
# quality assessment tool.
module Transrate

  Yell.new do |l|
    l.level = :info
    l.name = Object
    l.adapter STDOUT, level: [:debug, :info, :warn]
    l.adapter STDERR, level: [:error, :fatal]
  end
  Object.send :include, Yell::Loggable

end # Transrate
