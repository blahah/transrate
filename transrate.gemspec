# -*- encoding: utf-8 -*-
$:.push File.expand_path("../lib", __FILE__)
require File.expand_path('../lib/transrate/version', __FILE__)

Gem::Specification.new do |gem|
  gem.name          = 'transrate'
  gem.authors       = [ "Richard Smith-Unna", "Chris Boursnell", "Rob Patro", "Julian Hibberd", "Steven Kelly" ]
  gem.email         = "rds45@cam.ac.uk"
  gem.licenses      = ["MIT"]
  gem.homepage      = 'https://github.com/Blahah/transrate'
  gem.summary       = %q{ quality assessment of de-novo transcriptome assemblies }
  gem.description   = %q{ a library and command-line tool for quality assessment of de-novo transcriptome assemblies }
  gem.version       = Transrate::VERSION::STRING.dup

  gem.files = File.readlines('files.txt').map { |f| f.chomp }
  gem.executables = ["transrate"]
  gem.require_paths = %w( lib ext )
  gem.extensions  = ["ext/transrate/extconf.rb"]

  gem.add_dependency 'yell', '~> 2.0', '>= 2.0.4'
  gem.add_dependency 'trollop', '~> 2.0', '>= 2.0.0'
  gem.add_dependency 'bindeps', '~> 1.2', '>= 1.2.1'
  gem.add_dependency 'bio', '~> 1.4', '>= 1.4.3'
  gem.add_dependency 'crb-blast', '~> 0.6', '>= 0.6.4'
  gem.add_dependency 'colorize', '~> 0.7', '>= 0.7.7'

  gem.add_development_dependency 'test-unit', '~> 3.0'
  gem.add_development_dependency 'rake', '~> 10.3', '>= 10.3.2'
  gem.add_development_dependency 'rake-compiler', '~> 0.9', '>= 0.9.2'
  gem.add_development_dependency 'minitest', '~> 5'
  gem.add_development_dependency 'minitest-reporters', '~> 1', '>= 1.0.17'
  gem.add_development_dependency 'simplecov', '~> 0.8', '>= 0.8.2'
  gem.add_development_dependency 'shoulda', '~> 3.5', '>= 3.5.0'
  gem.add_development_dependency 'coveralls', '~> 0.7', '>= 0.7.2'
  gem.add_development_dependency 'bindeps', '~> 1.2', '>= 1.2.1'
end
