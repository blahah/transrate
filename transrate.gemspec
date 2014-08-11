# -*- encoding: utf-8 -*-
$:.push File.expand_path("../lib", __FILE__)
require File.expand_path('../lib/transrate/version', __FILE__)

Gem::Specification.new do |gem|
  gem.name          = 'transrate'
  gem.authors       = [ "Richard Smith-Unna", "Chris Boursnell" ]
  gem.email         = "rds45@cam.ac.uk"
  gem.licenses      = ["MIT"]
  gem.homepage      = 'https://github.com/Blahah/transrate'
  gem.summary       = %q{ quality assessment of de-novo transcriptome assemblies }
  gem.description   = %q{ a library and command-line tool for quality assessment of de-novo transcriptome assemblies }
  gem.version       = Transrate::VERSION::STRING.dup

  gem.files = `git ls-files`.split("\n")
  gem.executables = `git ls-files -- bin/*`.split("\n").map{ |f| File.basename(f) }
  gem.require_paths = %w( lib ext )
  gem.extensions  = ["ext/transrate/extconf.rb"]

  gem.add_dependency 'yell', '~> 2.0', '>= 2.0.4'
  gem.add_dependency 'trollop', '~> 2.0'
  gem.add_dependency 'bindeps', '~> 0.1', '>= 0.1.0'
  gem.add_dependency 'which', '~> 0.0', '>= 0.0.2'
  gem.add_dependency 'bio', '~> 1.4', '>= 1.4.3'
  gem.add_dependency 'crb-blast', '~> 0.4', '>= 0.4.1'
  gem.add_dependency 'bettersam', '~> 0.1', '>= 0.1.1'
  gem.add_dependency 'ruby-prof', '~> 0.15', '>= 0.15.1'

  gem.add_development_dependency 'rake', '~> 10.3', '>= 10.3.2'
  gem.add_development_dependency 'rake-compiler', '~> 0.9', '>= 0.9.2'
  gem.add_development_dependency 'turn', '~> 0.9', '>= 0.9.7'
  gem.add_development_dependency 'minitest', '~> 4', '>= 4.7.5'
  gem.add_development_dependency 'simplecov', '~> 0.8', '>= 0.8.2'
  gem.add_development_dependency 'shoulda-context', '~> 1.2', '>= 1.2.1'
  gem.add_development_dependency 'coveralls', '~> 0.7'
end
