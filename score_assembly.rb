#! /usr/bin/env ruby

require 'trollop'
require 'assembly'
require 'comparative_metrics'

opts = Trollop::options do
  version "v0.0.1a"
  banner <<-EOS

--------------------- k-descender ---------------------

run assembler (by default, velvet-oases) with descending
k-mers in a range (default, 90:25, step 6) for paired-end
reads

Richard Smith, April 2013

-------------------------------------------------------

EOS
  opt :assembly, "assembly file in FASTA format", :required => true, :type => String
  opt :reference, "reference proteome file in FASTA format", :required => true, :type => String
  opt :left, "left reads file in FASTQ format", :required => true, :type => String
  opt :right, "right reads file in FASTQ format", :required => true, :type => String
  opt :insertsize, "mean insert size",  :default => 200, :type => Integer
  opt :insertsd, "insert size standard deviation", :default => 50, :type => Integer
  opt :threads, "number of threads to use for parallel steps", :default => 8, :type => Integer
end

def pretty_print_hash hash, width
  hash.map{ |k, v| "#{k.to_s}#{" " * (width - (k.length + v.to_i.to_s.length))}#{v.to_i}" }.join("\n")
end

a = Assembly.new opts.assembly
r = Assembly.new opts.reference

puts "1. calculating basic stats..."
t0 = Time.now
basic_stats = a.basic_stats
puts "   ...done in #{Time.now - t0} seconds"

puts "\n2. calculating read diagnostics..."
t0 = Time.now
a.build_bowtie_index
mapped_reads = a.map_reads opts.left, opts.right, opts.insertsize, opts.insertsd
read_diagnostics = a.analyse_read_mappings mapped_reads, opts.insertsize, opts.insertsd
puts "   ...done in #{Time.now - t0} seconds"

puts "\n2. calculating comparative metrics..."
t0 = Time.now
comparative = ComparativeMetrics.new a, r
comparative_metrics = comparative.run
puts "   ...done in #{Time.now - t0} seconds"

puts "\n\n"
puts "|" *  40
puts "\n"
puts "Basic assembly metrics:"
puts "-" *  40
puts pretty_print_hash basic_stats, 40
puts "|" *  40
puts "\n"
puts "Read mapping metrics:"
puts "-" *  40
puts pretty_print_hash read_diagnostics, 40
puts "|" *  40
puts "\n"
puts "Comparative metrics:"
puts "-" *  40
puts pretty_print_hash comparative_metrics, 40
