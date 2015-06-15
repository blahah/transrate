require 'simplecov'
require 'coveralls'

SimpleCov.formatter = SimpleCov::Formatter::MultiFormatter[
  SimpleCov::Formatter::HTMLFormatter,
  Coveralls::SimpleCov::Formatter
]
SimpleCov.start

require 'minitest/autorun'
begin
  require 'turn/autorun'
  Turn.config.format = :pretty
  Turn.config.trace = 5
rescue LoadError
end
require 'shoulda/context'
require 'transrate'

# download large fastq files into test/data/.
path = "https://github.com/HibberdLab/transrate-test-files/raw/master"
test_files = []
test_files << File.join(File.dirname(__FILE__), "data", "sorghum_100.1.fastq")
test_files << File.join(File.dirname(__FILE__), "data", "sorghum_100.2.fastq")
test_files.each do |file|
  if !File.exist?(file)
    wget_cmd = "wget #{path}/#{File.basename(file)} --output-document #{file}"
    cmd = Transrate::Cmd.new(wget_cmd)
    cmd.run
  end
end
