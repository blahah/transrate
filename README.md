Transrate
----

Quality analysis and comparison of transcriptome assemblies.

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

**transrate** implements a variety of established and new metrics. 

note: this list will be expanded soon with detailed explanations and a guide to interpreting the results.

### Contig metrics

Contig metrics are measures based entirely on analysing the set of contigs themselves. At the moment, these are all related to the distribution of contig lengths and base content.

These are informative, but are only weakly useful for judging assembly quality. For most of these metrics, we don't know what the optimum is, although we can recognise extremely bad values. For example, an extremely small (\<1,000) or extremely large (\>500,000) number of contigs is biologically implausible for most organisms, and therefore suggests a problem with the assembly.

These metrics should therefore be used only as a quick, crude way of detecting major problems with the assembly.

| name          | explanation  |
| ------------- |:-------------|
| n_seqs | the number of contigs in the assembly |
| smallest | the size of the smallest contig |
| largest | the size of the largest contig |
| n_bases | the number of bases included in the assembly |
| mean_len | the mean length of the contigs |
| n > 1k | the number of contigs greater than 1,000 bases long |
| n > 10k | the number of contigs greater than 10,000 bases long |
| nX | the largest contig size at which at least X% of bases are contained in contigs at least this length |

### Read mapping metrics

Read mapping metrics are based on aligning a subset of the reads used in the assembly to the assembled contigs. These are a way of determining how well the assembly is supported by the original experimental evidence, and can be very useful overall quality metrics. However, they should not be considered alone, as read mapping metrics can be high for very fragmented assemblies.

These metrics can be useful for optimising your assembly. In particular, we want to maximise the proportion of the read pairs that map to the contigs successfully and in a biologically plausible way.

| name          | explanation  | optimum
| ------------- |:-------------| :----
| total | the total number and proportion of reads pairs mapping | theoretically 100%, but with erroneous and contaminating reads, often closer to 95%
| good | the number and proportion of read pairs mapping in a way indicative of good assembly | as above
| bad | the number and proportion of reads pairs mapping in a way indicative of bad assembly | 0%
| unexpressed contigs | the number and proportion of contigs that **are not** supported by read evidence | 0%
| expressed contigs | the number and proportion of contigs that **are** supported by read evidence | 100%

'Good' pairs are those aligned in a biologically plausible way, i.e.:

- where both members are aligned
- in the correct orientation
- either on the same contig or... 
- within a plausible distance of the ends of two separate contigs.

Conversely, 'bad' pairs are those where one of the conditions for being 'good' are not met.

Additionally, the software calculates whether there is any evidence in the read mappings that different contigs originate from the same transcript. These theoretical links are called bridges, and the number of bridges is shown in the **supported bridges** metric. A low count of supported bridges could be good or bad depending on the context. If you have a fragmented assembly, a large number of supported bridges means that scaffolding could greatly improve it. On the other hand, a large number of supported bridges in an otherwise seemingly good assembly could be indicative of misassemblies.

The list of supported bridges is output to a file, `supported_bridges.csv`, in case you want to make use of the information. At a later date, transrate will include the ability to scaffold the assembly using this and other information.

### Comparative metrics

| name          | explanation  | optimum
| ------------- |:-------------|:----
| reciprocal hits | the number of reciprocal best hits against the reference using ublast. A high score indicates that a large number of real transcripts have been assembled. | As high as possible. The theoretical maximum is the number of contigs (**n_seqs**). In practise, the maximum depends on the evolutionary divergence between the assembled species and the reference.
| contig hit rate | the proportion of contigs having a reciprocal best hit | As high as possible (see above)
| reference hit rate | the proportion of reference sequences having a reciprocal best hit | As high as possible (see above)
| ortholog hit ratio | the mean ratio of alignment length to reference sequence length. A low score on this metric indicates the assembly contains full-length transcripts. |  Close to 1
| collapse factor | the mean number of reference proteins mapping to each contig. A high score on this metric indicates the assembly contains chimeras. |  Dependent on the phylogenomic relationship between the organisms, e.g. whether a genome duplication has taken place.

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

`head -2000000 readfile.fastq`

FASTQ records are 4 lines long, so make sure you multiply the number of reads you want by 4, and be sure to run the same command on both the left and right read files.

### Example

```
transrate --assembly assembly.fasta \
	  --reference reference.fasta \
	  --left l.fq \
	  --right r.fq \
	  --threads 4
```

## Getting help

If you need help using transrate, please post to the [forum here](https://groups.google.com/forum/#!forum/transrate-users).

If you think you've found a bug, please post it to the [issues list](https://github.com/Blahah/transrate/issues).

## Requirements

### Ruby

First, you'll need Ruby v1.9.3 or greater installed. You can check with:

`ruby --version`

If you don't have Ruby installed, or you need a higher version, I recommend using [RVM](http://rvm.io/) as your Ruby Version Manager. To install RVM along with the latest Ruby, just run:

`\curl -L https://get.rvm.io | bash -s stable`

### Rubygems

Your Ruby installation *should* come with RubyGems, the package manager for Ruby. You can check with:

`gem --version`

If you don't have it installed, I recommend installing the latest version of Ruby and RubyGems using the RVM instructions above (in the Requirements:Ruby section).

### Usearch, Bowtie2 and eXpress

Usearch (http://drive5.com/usearch), Bowtie2 (https://sourceforge.net/projects/bowtie-bio/files/bowtie2) and eXpress (http://bio.math.berkeley.edu/eXpress/) must be installed and in your PATH. Additionally, the Usearch binary executable should be named `usearch`.
