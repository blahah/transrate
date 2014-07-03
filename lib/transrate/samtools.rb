require 'bio-samtools'

module Transrate

  class Samtools

    # Get the path to the samtools binary built when bio-samtools
    # was installed
    def self.path
      gem_path = Gem.loaded_specs['bio-samtools'].full_gem_path
      return File.join(gem_path, 'lib/bio/db/sam/external/samtools')
    end

    # Run a samtools command
    def self.run cmd
      runcmd = "#{Samtools.path} #{cmd}"
      runcmd
      `#{runcmd}`
    end

    # Convert a sam file to a bam file, returning the path to the bamfile
    def self.sam_to_bam samfile
      bamfile = File.basename(samfile, '.sam') + '.bam'
      bamfile
      Samtools.run "view -bS #{samfile} > #{bamfile}"
      File.expand_path bamfile
    end

    # Sort a bam file, returning the path to the sorted bamfile
    def self.sort_bam bamfile
      # the sort command behaves inconsistently with the other commands:
      # it takes an output prefix rather than a filename
      # and automatically adds the .bam extension
      sorted = File.basename(bamfile, '.bam') + '.sorted'
      Samtools.run "sort #{bamfile} #{sorted}"
      File.expand_path(sorted + '.bam')
    end

    # Index a bamfile, returning the path to the index
    def self.index_bam bamfile
      index = File.basename(bamfile, '.bam') + '.bai'
      Samtools.run "index #{bamfile} #{index}"
      File.expand_path index
    end

    def self.coverage(bam, contig)
      region = "#{contig.name}:1-#{contig.length}"
      result = Array.new(contig.length, 0)
      cmd = "mpileup"
      cmd += " -r #{region}" # region
      cmd += " -f #{bam.fasta}" # reference
      cmd += " -B" # don't calculate BAQ quality scores
      cmd += " -Q0" # include all reads ignoring quality
      cmd += " -I" # don't do genotype calculations
      cmd += " #{bam.bam}" # the bam file
      cov_idx = 3
      pos_idx = 1
      Samtools.run(cmd).split("\n").each do |line|
        cols = line.chomp.split("\t")
        cov = cols[cov_idx].to_i
        pos = cols[pos_idx].to_i
        result[pos - 1] = cov
      end
      result
    end

  end

end
