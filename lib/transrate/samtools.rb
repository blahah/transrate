module Transrate

  class Samtools

    class Bowtie2Error < StandardError; end

    require 'which'

    # Get the path to the samtools binary built when bio-samtools
    # was installed
    def self.path
      samtools_path = Which.which('samtools')
      if samtools_path.empty?
        raise SamtoolsError.new("could not find samtools in the path")
      end
      return samtools_path.kind_of?(String) ?
        samtools_path : samtools_path.first
    end

    # Run a samtools command
    def self.run cmd
      runcmd = Cmd.new "#{Samtools.path} #{cmd}"
      runcmd.run
      runcmd.stdout
    end

    # Convert a sam file to a bam file, returning the path to the bamfile
    def self.sam_to_bam samfile
      bamfile = File.basename(samfile, '.sam') + '.bam'
      bamfile = File.expand_path bamfile
      Samtools.run "view -bS #{File.expand_path samfile} > #{bamfile}"
      File.expand_path bamfile
    end

    # Sort a bam file, returning the path to the sorted bamfile
    def self.sort_bam bamfile
      # the sort command behaves inconsistently with the other commands:
      # it takes an output prefix rather than a filename
      # and automatically adds the .bam extension
      sorted = File.basename(bamfile, '.bam') + '.sorted'
      Samtools.run "sort #{File.expand_path bamfile} #{sorted}"
      File.expand_path(sorted + '.bam')
    end

    # Index a bamfile, returning the path to the index
    def self.index_bam bamfile
      index = File.basename(bamfile, '.bam') + '.bai'
      index = File.expand_path index
      Samtools.run "index #{File.expand_path bamfile} #{index}"
      index
    end

    # Convert a sam file to bam, sort and index the bam, returning
    # an array of paths to the bamfile, sorted bamfile and index respectively
    def self.sam_to_sorted_indexed_bam samfile
      bamfile = Samtools.sam_to_bam samfile
      sorted = Samtools.sort_bam bamfile
      index = Samtools.index_bam bamfile
      [bamfile, sorted, index]
    end

    # Calculate per-base coverage from a sorted, indexed bam file
    # return the path to the coverage file
    def self.coverage bam
      outfile = File.expand_path "#{File.basename(bam.fasta)}.coverage"
      cmd = "mpileup"
      cmd += " -f #{File.expand_path bam.fasta}" # reference
      cmd += " -B" # don't calculate BAQ quality scores
      cmd += " -Q0" # include all reads ignoring quality
      cmd += " -I" # don't do genotype calculations
      cmd += " #{File.expand_path bam.bam}" # the bam file
      cmd += " > #{outfile}"
      Samtools.run cmd
      outfile
    end

    # Calculate per-base coverage and mapQ score from a sorted, indexed
    # bam file. Return the path to the coverage file.
    def self.coverage_and_mapq(bam, fasta)
      outfile = File.expand_path "#{File.basename(fasta)}.bcf"
      cmd = "samtools mpileup"
      cmd << " -f #{File.expand_path fasta}" # reference
      cmd << " -B" # don't calculate BAQ quality scores
      cmd << " -q0" # include all multimapping reads
      cmd << " -Q0" # include all reads ignoring quality
      cmd << " -I" # don't do genotype calculations
      cmd << " -u" # output uncompressed bcf format
      cmd << " #{File.expand_path bam}" # the bam file
      cmd << " | bcftools view -cg - "
      cmd << " > #{outfile}"
      mpileup = Cmd.new cmd
      mpileup.run
      if !mpileup.status.success?
        raise RuntimeError.new("samtools and bcftools failed")
      end
      outfile
    end

  end

end
