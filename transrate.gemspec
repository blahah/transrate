# -*- encoding: utf-8 -*-
$:.push File.expand_path("../lib", __FILE__)
require File.expand_path('../lib/transrate/version', __FILE__)

Gem::Specification.new do |gem|
  gem.name          = 'transrate'
  gem.authors       = [ "Richard Smith" ]
  gem.email         = "rds45@cam.ac.uk"
  gem.homepage      = 'https://github.com/blahah/assemblotron'
  gem.summary       = %q{ quality assessment of de-novo transcriptome assemblies }
  gem.version       = Transrate::VERSION::STRING.dup

  gem.files = Dir['Rakefile', '{lib}/**/*', 'README*', 'LICENSE*']
  gem.require_paths = %w[ lib ]

  gem.add_dependency 'rake', '~> 10.1.0'
  gem.add_dependency 'trollop', '~> 2.0'
  gem.add_dependency 'which'
  gem.add_dependency 'bio'
  gem.add_dependency 'bettersam'

  gem.add_development_dependency 'turn'
  gem.add_development_dependency 'simplecov'
  gem.add_development_dependency 'shoulda-context'
  gem.add_development_dependency 'coveralls', '~> 0.6.7'
end