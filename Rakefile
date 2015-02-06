require 'rake/testtask'
require 'rake/extensiontask'

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
  t.name = :salmon
  t.libs << 'test'
  t.test_files = ['test/test_salmon.rb']
end


desc "Run tests"
task :default => :test

# PACKAGING

PACKAGE_NAME = "transrate"
VERSION = "1.0.0.beta2"
TRAVELING_RUBY_VERSION = "20141215-2.1.5"

desc "Package your app"
task :package => ['package:linux:x86_64', 'package:osx']

namespace :package do
  namespace :linux do
    desc "Package your app for Linux x86_64"
    task :x86_64 => "packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-linux-x86_64.tar.gz" do
      create_package("linux-x86_64")
    end
  end

  desc "Package your app for OS X"
  task :osx => "packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-osx.tar.gz" do
    create_package("osx")
  end
end

file "packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-linux-x86_64.tar.gz" do
  download_runtime("linux-x86_64")
end

file "packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-osx.tar.gz" do
  download_runtime("osx")
end

def create_package(target)
  package_dir = "packaging/#{PACKAGE_NAME}-#{VERSION}-#{target}"
  sh "rm -rf #{package_dir}"
  sh "mkdir -p #{package_dir}/lib/app"
  sh "cp -r lib #{package_dir}/lib/app/"
  sh "cp -r bin #{package_dir}/lib/app/"
  sh "mkdir #{package_dir}/lib/ruby"
  sh "tar -xzf packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-#{target}.tar.gz -C #{package_dir}/lib/ruby"
  sh "cp packaging/transrate #{package_dir}/transrate"
  if !ENV['DIR_ONLY']
    sh "tar -czf #{package_dir}.tar.gz #{package_dir}"
    sh "rm -rf #{package_dir}"
  end
end

def download_runtime(target)
  sh "mkdir -p packaging/packaging &&" +
  "cd packaging/packaging && curl -L -O --fail " +
  "http://d6r77u77i8pq3.cloudfront.net/releases/traveling-ruby-#{TRAVELING_RUBY_VERSION}-#{target}.tar.gz"
end
