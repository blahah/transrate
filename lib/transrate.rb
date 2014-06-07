# before the
require 'rbconfig'
RbConfig::CONFIG['CFLAGS'] = '-w'

require 'transrate/transrater'
require 'transrate/version'
require 'transrate/assembly'
require 'transrate/bowtie2'
require 'transrate/read_metrics'
require 'transrate/usearch'
require 'transrate/rb_hit'
require 'transrate/reciprocal_annotation'
require 'transrate/comparative_metrics'
require 'transrate/contig_metrics'
require 'transrate/metric'
require 'transrate/dimension_reduce'
require 'transrate/express'

# Transrate is a comprehensive transcriptome assembly
# quality assessment tool.
module Transrate

end # Transrate
