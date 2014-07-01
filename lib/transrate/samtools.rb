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
      bamfile
    end

    # Sort a bam file, returning the path to the sorted bamfile
    def self.sort_bam bamfile
      sorted = File.basename(bamfile, '.bam') + '.sorted'
      Samtools.run "sort #{bamfile} #{sorted}"
      sorted + '.bam'
    end

    # Index a bamfile, returning the path to the index
    def self.index_bam bamfile
      index = File.basename(bamfile, '.bam') + '.bai'
      Samtools.run "index #{bamfile} #{index}"
      index
    end

  end

end
