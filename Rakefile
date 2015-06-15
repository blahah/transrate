require 'rake/testtask'
require 'rake/extensiontask'
require 'bundler/setup'

Rake::ExtensionTask.new('transrate') do |ext|
  ext.lib_dir = "lib/transrate"
end

Rake::TestTask.new do |t|
  t.libs << 'test'
end

Rake::TestTask.new do |t|
  t.name = :recip
  t.libs << 'test'
  t.test_files = ['test/test_reciprocal.rb']
end

Rake::TestTask.new do |t|
  t.name = :comp
  t.libs << 'test'
  t.test_files = ['test/test_comp_metrics.rb']
end

Rake::TestTask.new do |t|
  t.name = :contig_metrics
  t.libs << 'test'
  t.test_files = ['test/test_contig_metrics.rb']
end

Rake::TestTask.new do |t|
  t.name = :read
  t.libs << 'test'
  t.test_files = ['test/test_read_metrics.rb']
end

Rake::TestTask.new do |t|
  t.name = :bowtie
  t.libs << 'test'
  t.test_files = ['test/test_bowtie.rb']
end

Rake::TestTask.new do |t|
  t.name = :rater
  t.libs << 'test'
  t.test_files = ['test/test_transrater.rb']
end

Rake::TestTask.new do |t|
  t.name = :bin
  t.libs << 'test'
  t.test_files = ['test/test_bin.rb']
end

Rake::TestTask.new do |t|
  t.name = :contig
  t.libs << 'test'
  t.test_files = ['test/test_contig.rb']
end

Rake::TestTask.new do |t|
  t.name = :assembly
  t.libs << 'test'
  t.test_files = ['test/test_assembly.rb']
end

Rake::TestTask.new do |t|
  t.name = :snap
  t.libs << 'test'
  t.test_files = ['test/test_snap.rb']
end

Rake::TestTask.new do |t|
  t.name = :salmon
  t.libs << 'test'
  t.test_files = ['test/test_salmon.rb']
end

Rake::TestTask.new do |t|
  t.name = :optimiser
  t.libs << 'test'
  t.test_files = ['test/test_optimiser.rb']
end



desc "Run tests"
task :default => :test

# PACKAGING

PACKAGE_NAME = "transrate"
VERSION = "1.0.0"
TRAVELING_RUBY_VERSION = "20150210-2.2.0"

desc "Package your app"
task :package => ['package:linux', 'package:osx']

namespace :package do
  desc "Package your app for Linux x86_64"
  task :linux => [:bundle_install, "packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-linux-x86_64.tar.gz"] do
    create_package("linux-x86_64")
  end

  desc "Package your app for OS X"
  task :osx => [:bundle_install, "packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-osx.tar.gz"] do
    create_package("osx")
  end
end

file "packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-linux-x86_64.tar.gz" do
  download_runtime("linux-x86_64")
end

file "packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-osx.tar.gz" do
  download_runtime("osx")
end

desc "Install gems to local directory"
task :bundle_install do
  if RUBY_VERSION !~ /^2\.2\./
    abort "You can only 'bundle install' using Ruby 2.2, because that's what Traveling Ruby uses."
  end
  Bundler.with_clean_env do
    sh "env BUNDLE_IGNORE_CONFIG=1 bundle install --path packaging/vendor --without development"
  end
  sh "rm -f packaging/vendor/*/*/cache/*"
end

def create_package(target)
  package_pref = "#{PACKAGE_NAME}-#{VERSION}-#{target}"
  package_dir = "packaging/#{package_pref}"
  sh "rm -rf #{package_dir}"
  sh "mkdir -p #{package_dir}/lib/app"
  # copy transrate gem to destination
  sh "cp -r lib bin deps ext files.txt #{package_dir}/lib/app/"
  # install travelling ruby
  sh "mkdir #{package_dir}/lib/app/ruby"
  sh "tar -xzf packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-#{target}.tar.gz -C #{package_dir}/lib/app/ruby"
  # install loading script for transrate
  sh "cp packaging/transrate #{package_dir}/transrate"
  # install bundled gem dependencies
  sh "cp -pR packaging/vendor #{package_dir}/lib/"
  sh "cd #{package_dir} && ../minify.sh"
  sh "cp -r #{package_dir}/lib/vendor/* #{package_dir}/lib/app/"
  sh "cp Gemfile Gemfile.lock transrate.gemspec #{package_dir}/lib/app/"
  sh "mkdir #{package_dir}/lib/app/.bundle"
  sh "cp packaging/bundler-config #{package_dir}/lib/app/.bundle/config"
  # free up some more space in the package dir
  sh "rm -rf #{package_dir}/lib/vendor"
  sh "rm -rf #{package_dir}/lib/app/ruby/*/gems/*/test"
  # install binary dependencies
  sh "mkdir -p packaging/bindeps/#{target}"
  sh "rm -rf packaging/bindeps/#{target}/*"
  sh "cp test/vagrant/#{target}/*.tar.gz packaging/bindeps/#{target}"
  sh "mkdir packaging/bindeps/#{target}/{bin,lib}"
  sh "cd packaging/bindeps/#{target} && " +
     "find . -maxdepth 1 -name '*.tar.gz' -exec tar -xzf '{}' \\; && " +
     "mv snap-aligner bam-read bin/"
  sh "cp -r packaging/bindeps/#{target}/{bin,lib} #{package_dir}/"
  # install c extension
  sh "cp test/vagrant/#{target}/{transrate,libruby}.* #{package_dir}/lib/"
  # create package
  if !ENV['DIR_ONLY']
    sh "cd packaging && tar -czf #{package_pref}.tar.gz #{package_pref}"
    sh "rm -rf #{package_dir}"
  end
  # cleanup
  sh "rm -rf packaging/vendor packaging/bindeps .bundle"
end

def download_runtime(target)
  sh "mkdir -p packaging/packaging &&" +
  "cd packaging/packaging && curl -L -O --fail " +
  "http://d6r77u77i8pq3.cloudfront.net/releases/traveling-ruby-#{TRAVELING_RUBY_VERSION}-#{target}.tar.gz"
end
