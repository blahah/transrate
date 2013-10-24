Transrate
----

Quality analysis and comparison of transcriptome assemblies.

## Contents

1. [Development status](https://github.com/Blahah/transrate#development-status)
2. [Transcriptome assembly quality metrics](https://github.com/Blahah/transrate#transcriptome-assembly-quality-metrics)
3. [Installation](https://github.com/Blahah/transrate#installation)
4. [Usage](https://github.com/Blahah/transrate#usage)
   - [example](https://github.com/Blahah/transrate#example)
5. [Requirements](https://github.com/Blahah/transrate#requirements)
   1. [Ruby](https://github.com/Blahah/transrate#ruby)
   2. [RubyGems](https://github.com/Blahah/transrate#rubygems)
   3. [USEARCH, Bowtie 2, and eXpress](https://github.com/Blahah/transrate#usearch-bowtie2-and-express)
6. [Getting help](https://github.com/Blahah/transrate#getting-help)

## Development status

This software is in early development. Users should be aware that until the first release is made, features may change faster than the documentation is updated. Nevertheless, we welcome bug reports.

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

## Transcriptome assembly quality metrics

**transrate** implements a variety of established and new metrics. They are explained in detail [on the wiki](https://github.com/Blahah/transrate/wiki/Transcriptome-assembly-quality-metrics).

## Installation

Assuming all the requirements are met (see below), you can install transrate very easily. Just run at the terminal:

`gem install transrate`

If you're new to linux/unix, there's a detailed tutorial for installing transrate with all the dependencies [on my blog](http://blahah.net/bioinformatics/2013/10/19/installing-transrate/).

## Usage

`transrate --help` will give you...

```
Transrate v0.0.10 by Richard Smith <rds45@cam.ac.uk>

DESCRIPTION:
Analyse a de-novo transcriptome
assembly using three kinds of metrics:

1. contig-based
2. read-mapping
3. reference-based

Please make sure USEARCH, bowtie 2 and eXpress are installed
and in the PATH.

Bug reports and feature requests at:
http://github.com/blahah/transrate

USAGE:
transrate <options>

OPTIONS:
    --assembly, -a <s>:   assembly file in FASTA format
   --reference, -r <s>:   reference proteome file in FASTA format
        --left, -l <s>:   left reads file in FASTQ format
       --right, -i <s>:   right reads file in FASTQ format
  --insertsize, -n <i>:   mean insert size (default: 200)
    --insertsd, -s <i>:   insert size standard deviation (default: 50)
     --threads, -t <i>:   number of threads to use (default: 8)
         --version, -v:   Print version and exit
            --help, -h:   Show this message
```

If you don't include --left and --right read files, the read-mapping based analysis will be skipped. I recommend that you don't align all your reads - just a subset of 500,000 will give you a very good idea of the quality. You can get a subset by running (on a linux system):

`head -2000000 left.fastq > left_500k.fastq`
`head -2000000 right.fastq > right_500k.fastq`

FASTQ records are 4 lines long, so make sure you multiply the number of reads you want by 4, and be sure to run the same command on both the left and right read files.

### Example

```
transrate --assembly assembly.fasta \
	  --reference reference.fasta \
	  --left l.fq \
	  --right r.fq \
	  --threads 4
```

## Requirements

### Ruby

First, you'll need Ruby v1.9.3 or greater installed. You can check with:

`ruby --version`

If you don't have Ruby installed, or you need a higher version, I recommend using [RVM](http://rvm.io/) as your Ruby Version Manager. To install RVM along with the latest Ruby, just run:

`\curl -L https://get.rvm.io | bash -s stable`

### Rubygems

Your Ruby installation *should* come with RubyGems, the package manager for Ruby. You can check with:

`gem --version`

If you don't have it installed, I recommend installing the latest version of Ruby and RubyGems using the RVM instructions above (in the [Requirements:Ruby](https://github.com/Blahah/transrate#ruby) section).

### Usearch, Bowtie2 and eXpress

Usearch (http://drive5.com/usearch), Bowtie2 (https://sourceforge.net/projects/bowtie-bio/files/bowtie2) and eXpress (http://bio.math.berkeley.edu/eXpress/) must be installed and in your PATH. Additionally, the Usearch binary executable should be named `usearch`.

## Getting help

If you need help using transrate, please post to the [forum here](https://groups.google.com/forum/#!forum/transrate-users).

If you think you've found a bug, please post it to the [issues list](https://github.com/Blahah/transrate/issues).
