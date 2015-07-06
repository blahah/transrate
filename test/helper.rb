require 'simplecov'
require 'coveralls'

SimpleCov.formatter = SimpleCov::Formatter::MultiFormatter[
  SimpleCov::Formatter::HTMLFormatter,
  Coveralls::SimpleCov::Formatter
]
SimpleCov.start

require "stringio"

def capture_stderr
  real_stderr, $stderr = $stderr, StringIO.new
  yield
  $stderr.string
ensure
  $stderr = real_stderr
end

def capture_stdout
  real_stdout, $stdout = $stdout, StringIO.new
  yield
  $stdout.string
ensure
  $stdout = real_stdout
end

# use within an at_exit block, exploits the capture of last exception in
# the global var $!
def last_exit_successful?
  $!.nil? || $!.is_a?(SystemExit) && $!.success?
end

def sorghum_data
  assembly = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.fa')
  reference = File.join(File.dirname(__FILE__), 'data', 'Os.protein.2.fa')
  left = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.1.fastq')
  right = File.join(File.dirname(__FILE__), 'data', 'sorghum_100.2.fastq')
  [assembly, reference, left, right]
end

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
