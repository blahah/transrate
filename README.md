Transrate
----

Quality analysis and comparison of transcriptome assemblies.

[![Gem Version](https://badge.fury.io/rb/transrate.png)][gem]
[![Build Status](https://secure.travis-ci.org/Blahah/transrate.png?branch=master)][travis]
[![Dependency Status](https://gemnasium.com/Blahah/transrate.png?travis)][gemnasium]
[![Code Climate](https://codeclimate.com/github/Blahah/transrate.png)][codeclimate]
[![Coverage Status](https://coveralls.io/repos/Blahah/transrate/badge.png?branch=master)][coveralls]

[gem]: https://badge.fury.io/rb/transrate
[travis]: https://travis-ci.org/Blahah/transrate
[gemnasium]: https://gemnasium.com/Blahah/transrate
[codeclimate]: https://codeclimate.com/github/Blahah/transrate
[coveralls]: https://coveralls.io/r/Blahah/transrate

## Contents

1. [Development status](https://github.com/Blahah/transrate#development-status)
2. [Transcriptome assembly quality metrics](https://github.com/Blahah/transrate#transcriptome-assembly-quality-metrics)
3. [Installation](https://github.com/Blahah/transrate#installation)
4. [Usage](https://github.com/Blahah/transrate#usage)
    - [Command line](https://github.com/Blahah/transrate#command-line)
        - [example](https://github.com/Blahah/transrate#example)
    - [As a library](https://github.com/Blahah/transrate#as-a-library)
5. [Requirements](https://github.com/Blahah/transrate#requirements)
    - [Ruby](https://github.com/Blahah/transrate#ruby)
    - [RubyGems](https://github.com/Blahah/transrate#rubygems)
    - [Blast+, Bowtie 2](https://github.com/Blahah/transrate#blast+-and-bowtie2)
6. [Getting help](https://github.com/Blahah/transrate#getting-help)

## Development status

This software is being actively developed. Please be aware that they may be bugs. If you find any, please report them on the [issue tracker](https://github.com/Blahah/transrate/issues).

## Transcriptome assembly quality metrics

**transrate** implements a variety of established and new metrics. They are explained in detail [on the wiki](https://github.com/Blahah/transrate/wiki/Transcriptome-assembly-quality-metrics).

## Installation

Assuming you've got a recent version of Ruby installed (see below), you can install transrate very easily. Just run at the terminal:

`gem install transrate`

Next all the software Transrate depends on needs to be installed. Luckily, transrate is clever enough to do this itself. Simply run

`transrate --install-deps`

Transrate will check whether its dependencies are installed, and if not will download and install them for you.

## Usage

### Command line

`transrate --help` will display basic usage instructions.

```
  Transrate v0.2.0 by Richard Smith-Unna <rds45@cam.ac.uk>

  DESCRIPTION:
  Analyse a de-novo transcriptome
  assembly using three kinds of metrics:

  1. contig-based
  2. read-mapping (if --left and --right are provided)
  3. reference-based (if --reference is provided)

  Bug reports and feature requests at:
  http://github.com/blahah/transrate

  USAGE:
  transrate <options>

  EXAMPLES:
  transrate --assembly contigs.fa --reference Athaliana_protein.fa --threads 8

  OPTIONS:
    --assembly, -a <s>:   assembly file(s) in FASTA format, comma-separated
   --reference, -r <s>:   reference proteome file in FASTA format
        --left, -l <s>:   left reads file in FASTQ format
       --right, -i <s>:   right reads file in FASTQ format
  --insertsize, -n <i>:   mean insert size (default: 200)
    --insertsd, -s <i>:   insert size standard deviation (default: 50)
     --threads, -t <i>:   number of threads to use (default: 8)
     --outfile, -o <s>:   filename to use for CSV output (default: transate_results.csv)
    --loglevel, -g <s>:   the amount of information to print. one of [error, info, warn, debug] (default: info)
    --install-deps, -d:   install any missing dependencies
         --profile, -p:   debug option: profile the code as it runs
         --version, -v:   Print version and exit
            --help, -h:   Show this message
```

See the [getting started guide] on the website for more instructions, and see the [command-line options] part of the manual for details.

#### Example

```
transrate --assembly assembly.fasta \
	  --reference reference.fasta \
	  --left l.fq \
	  --right r.fq \
	  --threads 4
```

### As a library

```ruby
require 'transrate'

assembly = Transrate::Assembly.new(File.expand_path('assembly.fasta'))
reference = Transrate::Assembly.new(File.expand_path('reference.fasta'))

t = Transrate::Transrater.new(assembly, reference)

left = File.expand_path('left.fq')
right = File.expand_path('right.fq')

puts t.all_metrics(left, right)
puts t.assembly_score
```

## Requirements

### Ruby

First, you'll need Ruby v2.0.0 or greater installed. You can check with:

`ruby --version`

If you don't have Ruby installed, or you need a higher version, I recommend using [RVM](http://rvm.io/) as your Ruby Version Manager. To install RVM along with the latest Ruby, just run:

`\curl -L https://get.rvm.io | bash -s stable --ruby`

### Rubygems

Your Ruby installation *should* come with RubyGems, the package manager for Ruby. You can check with:

`gem --version`

If you don't have it installed, I recommend installing the latest version of Ruby and RubyGems using the RVM instructions above (in the [Requirements:Ruby](https://github.com/Blahah/transrate#ruby) section).

## Other software

Transrate uses a variety of other software including NCBI BLAST+, Bowtie2 and Samtools. But you don't need to worry about installing those yourself - Transrate does that automatically using the [bindeps](https://github.com/Blahah/bindeps) gem.

## Getting help

If you need help using transrate, please post to the [forum here](https://groups.google.com/forum/#!forum/transrate-users).

If you think you've found a bug, please post it to the [issues list](https://github.com/Blahah/transrate/issues).
