# before the
require 'rbconfig'
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
require 'transrate/express'
require 'transrate/samtools'

# Transrate is a comprehensive transcriptome assembly
# quality assessment tool.
module Transrate

end # Transrate
