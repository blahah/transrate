Transrate
----

Quality analyis and comparison of transcriptome assemblies.


## Installation

You can install transrate very easily. Just run at the terminal:

`gem install transrate`

If that doesn't work, check the requirements below...

## Usage

`transrate --help` will give you...

```
Transrate v0.0.1a by Richard Smith <rds45@cam.ac.uk>

DESCRIPTION:
Analyse a de-novo transcriptome
assembly using three kinds of metrics:

1. contig-based
2. read-mapping
3. reference-based

Please make sure USEARCH and bowtie2 are both installed
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

If you don't have it installed, I recommend installing the latest version of Ruby and RubyGems using the RVM instructions above (in the Requirements:Ruby section.

## Development status

This software is in very early development. Nevertheless, we welcome bug reports.