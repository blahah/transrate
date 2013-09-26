#!/usr/bin/env ruby
# take a transcriptome as input, and output replicas with different
# types of errors added at different rates

# TODO:
# - allow variable insertion/deletion size
# - include non-self chimeras
# - collect more error types

require 'bio'
require 'rubystats'

class Bio::Sequence::NA
  
  # return a chimera
  # if other_seq=nil, chimera is created by
  # binding the sequence to itself.
  def chimera(other_seq=nil)
    other_seq ||= self.randomflip
    self + other_seq
  end

  # return sequence randomly flipped
  # in orientation and complementarity
  def randomflip
    seq = self
    r = rand
    if r < 0.25
      seq = seq.reverse
    elsif r < 0.5
      seq = seq.reverse_complement
    elsif r < 0.75
      seq = seq.forward_complement
    end
    Bio::Sequence::NA.new(seq)
  end

  # return two sequences
  # created by fragmenting the sequence at
  # a point drawn from a normal distribution
  # centered around the middle of the sequence,
  # and with a standard deviation of 1/5 the
  # sequence length
  def fragment
    pos = self.rnorm_position(self.length)
    seq = self.to_s
    s1 = seq[0..pos]
    s2 = seq[pos + 1..-1]
    [Bio::Sequence::NA.new(s1), Bio::Sequence::NA.new(s2)]
  end

  # given a sequence, return a sequence with
  # a random insertion of size at a point drawn
  # froma  normal distribution centered around
  # the middle of the sequence, and with a SD
  # of 1/5 the sequence length
  def insertion size=5
    pos = self.rnorm_position(self.length)
    newseq = self[0..size].randomize
    self[0..size] + newseq + self[size + 1..-1]
  end

  # given a sequence, return a sequence with
  # a deletion of size centered around a point
  # drawn with rnorm_position
  def deletion size=5
    pos = self.rnorm_position(self.length)
    self[0..pos - (size / 2)] + self[pos + (size / 2) + 1..-1]
  end

  # given a sequence, return two sequences created
  # by simulated bubbling of an assembly graph at
  # a point drawn with rnorm_position
  def bubble
    pos = self.rnorm_position(self.length)
    letter = case self[pos]
             when "a"
               "t"
             when "t"
               "a"
             when "g"
               "c"
             when "c"
               "g"
             end
    seq = self.to_s
    seq[pos] = letter
    [self, Bio::Sequence::NA.new(seq)]
  end

  # given a length, return a number drawn
  # from a normal distrubition with mean
  # length/2 and SD length/5
  def rnorm_position(length)
    p = Rubystats::NormalDistribution.new(length.to_f / 2, length.to_f / 6).rng
    p = length / 2 if p >= length
    p = length / 2 if p <= 0
    p
  end

end # Bio::Sequence::NA


require 'trollop'
opts = Trollop::options do
  opt :file, "transcriptome FASTA file to simulate errors from", :type => :string
end

$sequences = []
Bio::FastaFormat.open(opts.file).each do |entry|
  $sequences << entry.naseq
end

def simulate_txome(filename, rate, method)
  i = 0
  File.open(filename, 'w') do |f|
    $sequences.each do |seq|
      r = rand
      seq = seq.send(method) if r < rate
      if seq.is_a? Array
        seq.each_with_index do |s, ii|
          f.puts s.to_fasta("#{i}_#{ii}", 60)
        end
      else
        f.puts seq.to_fasta("#{i}", 60)
      end
    end
    i += 1
  end
end

rates = [0.01, 0.05, 0.10, 0.25, 0.5]
modifiers = [:chimera, :fragment, :insertion, :deletion, :bubble]
rates.each do |rate|
  modifiers.each do |mod|
    puts "simulating transcriptome with #{mod.to_s} at rate #{rate}"
    simulate_txome("#{mod.to_s}_#{rate}.fa", rate, mod)
  end
end
