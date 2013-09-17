module Transrate

  class RBHit

    # Fields: query id, subject id, % identity, alignment length, mismatches,
    # gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    attr_accessor :query, :target, :id, :alnlen, :mismatches
    attr_accessor :gaps, :qstart, :qend, :tstart, :tend, :evalue
    attr_accessor :bitscore, :tcov
    
    def initialize(list)
      @query      = list[0].scan(/[^|]+/).first.split.first # extract only identifier
      @target     = list[1].scan(/[^|]+/).first.split.first
      @id         = list[2]
      @alnlen     = list[3]
      @mismatches = list[4]
      @gaps       = list[5]
      @qstart     = list[6]
      @qend       = list[7]
      @tstart     = list[8]
      @tend       = list[9]
      @evalue     = list[10]
      @bitscore   = list[11]
      @tcov       = list[12]
    end

    def to_s
      @query + " => " + @target
    end

  end # RBHit

end # Transrate