#!/usr/bin/env ruby
# This is a short test script to demonstrate that SAM files
# generated with the current bowtie2 parameters
# do not repeat reads, even those with multiple alignments.
# The sam file includes the exact command used to run bowtie
test={}
samfile = File.join(File.dirname(__FILE__), "test.sam")
File.open(samfile,'r').each do |line|
  test[line.split[0]] ? test[line.split[0]] += 1 : test[line.split[0]] = 1
end
test.each do |key,value|
  raise "more than one instance of a read in a pair" if value > 2
end

puts "All Finished! No repeats detected"
