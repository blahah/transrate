# before the
require 'rbconfig'
RbConfig::CONFIG['CFLAGS'] = ''

require 'transrate/log'
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

  def self.log
    Log.instance
  end

end # Transrate
