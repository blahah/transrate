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
  t.name = :refaln
  t.libs << 'test'
  t.test_files = ['test/test_reference_alignment.rb']
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

desc "Run tests"
task :default => :test
