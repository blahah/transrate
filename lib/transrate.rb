require 'transrate/transrater'
require 'transrate/version'
require 'transrate/assembly'
require 'transrate/bowtie2'
require 'transrate/read_metrics'
require 'transrate/usearch'
require 'transrate/rb_hit'
require 'transrate/reciprocal_annotation'
require 'transrate/comparative_metrics'
require 'transrate/metric'
require 'transrate/dimension_reduce'
require 'transrate/express'
require 'transrate/log'

module Transrate

  # one logger instance is stored as a module ivar
  # to avoid creating lots of loggers
  class << self
    attr_accessor :log
  end # log accessor

end # Transrate
