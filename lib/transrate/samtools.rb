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
      `#{runcmd}`
    end

  end

end
