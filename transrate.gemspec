# -*- encoding: utf-8 -*-
$:.push File.expand_path("../lib", __FILE__)
require File.expand_path('../lib/transrate/version', __FILE__)

Gem::Specification.new do |gem|
  gem.name          = 'transrate'
  gem.authors       = [ "Richard Smith" ]
  gem.email         = "rds45@cam.ac.uk"
  gem.licenses      = ["MIT"]
  gem.homepage      = 'https://github.com/Blahah/transrate'
  gem.summary       = %q{ quality assessment of de-novo transcriptome assemblies }
  gem.description   = %q{ a library and command-line tool for quality assessment of de-novo transcriptome assemblies }
  gem.version       = Transrate::VERSION::STRING.dup

  gem.files = `git ls-files`.split("\n")
  gem.executables = `git ls-files -- bin/*`.split("\n").map{ |f| File.basename(f) }
  gem.require_paths = %w( lib )

  gem.add_dependency 'rake'
  gem.add_dependency 'trollop'
  gem.add_dependency 'which'
  gem.add_dependency 'bio'
  gem.add_dependency 'bio-samtools'
  gem.add_dependency 'RubyInline', '~> 3.12'
  gem.add_dependency 'crb-blast'

  gem.add_development_dependency 'test-unit'
  gem.add_development_dependency 'minitest'
  gem.add_development_dependency 'turn'
  gem.add_development_dependency 'simplecov'
  gem.add_development_dependency 'shoulda-context'
  gem.add_development_dependency 'coveralls', '>= 0.6.7'
end
