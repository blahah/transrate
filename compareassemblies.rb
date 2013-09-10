#!/usr/bin/env ruby

require 'rubygems'
require 'trollop'
require 'bio'
require 'zlib'
require 'csv'
require 'jruby_threach'

# TODO: add thread options
# TODO: overhaul the script:
#   - keep this file as a controller which runs the components and reports progress
#   - split out other functions into individual scripts (e.g. usearcher, reference assembly metrics, merge all assembly metrics, )

# options
opts = Trollop::options do
  version "v0.0.1a"
  banner <<-EOS

----------------- compare_assemblies ------------------

compare_assemblies analyses transcriptome assemblies in 
comparison to a provided reference. Outputs a report of
meaningful metrics to allow a deep relative quality 
assessment between assemblies.

--------------------------------------------------------

EOS
  opt :ref, "Reference assembly file name (required)", :required => true, :type => String
  opt :refexpression, "Reference expression quantification (RSEM) file name (required)", :required => true, :type => String
  opt :inputs, "List of input assembly file names [e.g. 1.fa,2.fa,5.fa,10.fa](required)", :required => true, :type => String
  opt :labels, "List of labels to assign to input assemblies [e.g. 1mil,2mil,5mil,10mil](required)", :required => true, :type => String
  opt :sid, "usearch sequence identity threshold as float (defaults to 0.95)", :default => 0.95, :type => Float
  opt :onlyusearch, "Just perform the usearch, then exit"
  opt :skipusearch, "Skip the usearch, just analyse the provided output"
  opt :debug, "Be extremely verbose for debugging purposes"

end

TESTASSEMBLYNAMES = opts.inputs.split(",")
TESTASSEMBLYLABELS = opts.labels.split(",")
REFERENCEASSEMBLY = opts.ref
REFERENCEEXPRESSIONFILE = opts.refexpression
SEARCHSIMILARITY = opts.sid
DEBUG = opts.debug

# helper functions
def execOrError(cmd)
  # execute a system command, catching the exit status
  # and reporting any errors to the user
  retvalue = `#{cmd}`
  # exit status object is automatically stored as $?
  # exit code of 0 is good, everything else bad
  if $?.exitstatus != 0
    puts "command '#{cmd}' failed with error: #{retvalue}"
  end
  return retvalue
end

def usearchQuery(assembly, precoverage, scope='global')
  return "~/apps/usearch -usearch_#{scope} #{assembly} -db reference.udb --strand both -threads 8" \
  " -id #{SEARCHSIMILARITY} -userout #{assembly}.#{scope} #{precoverage} -fastapairs aligned.#{assembly}" \
  " -matched matched.#{assembly}.#{scope} -dbmatched matched.#{assembly}.#{REFERENCEASSEMBLY}.#{scope}"
end

class Bio::Sequence
  def complexity
    sequence = self.to_s
    compressed = Zlib::Deflate.deflate sequence
    return compressed.length.to_f/sequence.length
  end

  def entropy
    # from https://gist.github.com/747242
    len = self.to_s.chars.count.to_f
    log2 = Math.log(2)
    
    counts = word.chars.inject({}) do |h,c|
      h[c] = (h[c] || 0) + 1
      h
    end
   
    counts.inject(0) do |entropy, pair|
      frequency = (pair[1] / len)
      entropy = (frequency * (Math.log(frequency) / log2))
    end.abs

    return entropy
  end
end

def writeTSV(filename, array)
  CSV.open(filename, 'w', {:col_sep => "\t"}) do |tsv|
    array.each do |line|
      tsv << line
    end
  end
end

=begin
test gene-level coverage
usearch_global each assembly against reference
this gives us the blast output, plus the fasta files
of genes matched in both the reference and the assembly
=end

unless opts.skipusearch
  # make the reference database
  if !File.exists? "reference.udb"
    puts 'building reference usearch database...'
    execOrError "~/apps/usearch -makeudb_usearch #{REFERENCEASSEMBLY} -output reference.udb"
  end

  # search each assembly into the reference
  precoverage = '-userfields query+target+id+tcov'
  TESTASSEMBLYNAMES.threach(3) do |assembly|
    if !File.exists? "matched.#{assembly}.#{REFERENCEASSEMBLY}.global" or
        !File.exists? "matched.#{assembly}.global"
      puts '-' * 40
      puts "searching #{assembly} against reference..."
      query = usearchQuery(assembly, precoverage)
      puts query
      execOrError query
    end
  end
end


puts 'usearching done!'

if opts.onlyusearch
  exit
end

=begin
now for each assembly we have:
blast-format output: assemblyname.b6 
contigs found in reference: matched.assemblyname.global
transcripts found in assembly: matched.assemblyname.referencename.global
=end

puts '-' * 40
puts 'preparing to calculate per-macromolecule metrics...'
$countcommand = "grep -c '^>' "
$refcount = `#{$countcommand} #{REFERENCEASSEMBLY}`.to_f
puts "refcount = #{$refcount}"
$totalcounts = [["reference", $refcount]]
$globalvsreferenceproportions = [] # [name, count] for each assembly
$globalvsassemblyproportions = [] # [name, count] for each assembly
  
# macro coverage for each assembly
TESTASSEMBLYNAMES.each_with_index do |assembly, index|
  label = TESTASSEMBLYLABELS[index]
  assemblycount = `#{$countcommand} #{assembly}`.to_f
  globaltranscriptsmatched = `#{$countcommand} matched.#{assembly}.#{REFERENCEASSEMBLY}.global`.to_f
  globalcontigsmatched = `#{$countcommand} matched.#{assembly}.global`.to_f
  $totalcounts << [label, assemblycount]
  $globalvsreferenceproportions << [label, globaltranscriptsmatched / $refcount]
  $globalvsassemblyproportions << [label, globalcontigsmatched / assemblycount]
end

puts 'calculations done!'
puts '-' * 40
puts 'generating data output...'
output = 'analysis'

# write out total counts bar chart data
writeTSV(output+'/totalcounts.tsv', $totalcounts)
# write out total reference coverage bar chart data
writeTSV(output+'/globalreferencecoverage.tsv', $globalvsreferenceproportions)
# write out total assembly coverage bar chart data
writeTSV(output+'/globalassemblycoverage.tsv', $globalvsassemblyproportions)

$transcripts = {}

# for testing, check if we've already got a saved transcript hash
unless File.exist?('transcripts.final.hash')
  unless File.exist?('transcripts.hash')
    # REFERENCE
    puts 'parsing reference assembly'
    # synchronously step through the reference and expression
    afile = Bio::FastaFormat.open("#{REFERENCEASSEMBLY}") # assembly
    efile = CSV.open("#{REFERENCEEXPRESSIONFILE}", 'r', {:col_sep => ' '}) # expression
    efile.readline # discard the header
    expr = efile.readline # current line of expression
    afile.each do |fasta|
      # check the ids match 
      fpkm = 0.0
      if fasta.definition == expr[0]
        fpkm = expr[6].to_f
        expr = efile.readline
      end
      s = Bio::Sequence.auto(fasta.seq)
      # add the metrics for the reference transcript
      transcriptdata = {
        :l => s.length,
        :g => s.gc_content * 100,
        :c => s.complexity,
        :e => fpkm
      }
      # store it in the hash
      $transcripts[fasta.definition] = transcriptdata
    end
    puts 'reference assembly measurements complete'

    # for testing purposes, dump transcripts to file
    File.open('transcripts.hash', 'w') {|f| f.write(Marshal.dump($transcripts)) }
  else
    # load saved transcript hash
    puts 'loading saved reference transcript hash'
    $transcripts = Marshal.load(File.read('transcripts.hash'))
  end
  puts 'saving raw transcript hash...'
  File.open('transcripts.final.hash', 'w') {|f| f.write(Marshal.dump($transcripts)) }
  puts '...done'
else
  puts 'loading saved final transcript hash'
  $transcripts = Marshal.load(File.read('transcripts.final.hash'))
end

puts 'outputting data table'

# column headings so we can label outputs
$header = ['Transcript ID','Length','GC%','Complexity','FPKM']
# construct a data table
CSV.open("metricstable.csv", "w") do |table|
  table << $header
  $transcripts.each do |transcript, data|
    line = [transcript]
    line << data[:l]
    line << data[:g]  
    line << data[:c]
    line << data[:e] # expression
    table << line
  end
end

puts 'outputting sample data'

CSV.open("samples.csv", "w") do |table|
  table << ['filename', 'label', 'searchresults', 'alignment']
  (0..TESTASSEMBLYNAMES.length-1).each do |i|
    line = [TESTASSEMBLYNAMES[i]]
    line << TESTASSEMBLYLABELS[i]
    line << "#{line[0]}.global"
    line << "aligned.#{line[0]}"
    table << line
  end
end

puts 'done!'
