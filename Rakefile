require 'rake/testtask'

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

desc "Run tests"
task :default => :test