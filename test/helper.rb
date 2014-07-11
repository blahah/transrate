require 'simplecov'
require 'coveralls'

SimpleCov.formatter = SimpleCov::Formatter::MultiFormatter[
  SimpleCov::Formatter::HTMLFormatter,
  Coveralls::SimpleCov::Formatter
]
SimpleCov.start

require 'minitest/autorun'
begin; require 'turn/autorun'; rescue LoadError; end
require 'shoulda/context'
require 'transrate'

Turn.config.format = :pretty
Turn.config.trace = 5

# fake CRBBlast class
class CRBHelper

  attr_accessor :target_is_prot, :hash
  def initialize t
    @target_is_prot = t
  end

  def reciprocals
    return @hash
  end

end

# rake Hit class
class HitHelper

  attr_accessor :query, :target, :qstart, :qend, :tstart, :tend, :qlen, :tlen
  def initialize query, target, qstart, qend, tstart, tend, qlen, tlen
    @query = query
    @target = target
    @qstart = qstart
    @tstart = tstart
    @tend = tend
    @qend = qend
    @qlen = qlen
    @tlen = tlen
  end

end
