module Transrate

class Cmdline

  require 'trollop'
  require 'csv'
  require 'bindeps'
  require 'colorize'
  require 'pathname'

  def initialize args
    @opts = parse_arguments args
    if @opts.examples
      print_examples
    end
    @report_width = 35
    @outfile = "assemblies.csv"
    check_arguments
  end

  def run
    results = []

    assemblies = @opts.assembly.split(',')
    result_paths = assembly_result_paths assemblies

    r = @opts.reference ? Assembly.new(File.expand_path @opts.reference) : nil

    @opts.output = File.expand_path @opts.output
    FileUtils.mkdir_p @opts.output

    Dir.chdir @opts.output do
      if @opts.merge_assemblies
        assemblies = concatenate_assemblies assemblies
      end

      assemblies.zip(result_paths) do |assembly, result_path|
        results << analyse_assembly(assembly, r, result_path)
      end

      write_assembly_csv results
    end

  end

  def parse_arguments args
    Trollop::with_standard_exception_handling argument_parser do
      if args.empty? || args.include?("-h") || args.include?("--help")
        transrate_banner
        raise Trollop::HelpNeeded
      end

      argument_parser.parse args
    end
  end

  def argument_parser
    cmdline = self
    Trollop::Parser.new do
      version Transrate::VERSION::STRING.dup
      banner cmdline.help_message
      opt :assembly, "Assembly file(s) in FASTA format, comma-separated",
          :type => String
      opt :left, "Left reads file(s) in FASTQ format, comma-separated",
          :type => String
      opt :right, "Right reads file(s) in FASTQ format, comma-separated",
          :type => String
      opt :reference,
          "Reference proteome or transcriptome file in FASTA format",
          :type => String
      opt :threads, "Number of threads to use",
          :default => 8,
          :type => Integer
      opt :merge_assemblies,
          "Merge best contigs from multiple assemblies into file",
          :type => String
      opt :output, "Directory where results are output (will be created)",
          :default => 'transrate_results'
      opt :loglevel,
          "Log level. One of [error, info, warn, debug]",
          :default => 'info'
      opt :install_deps,
          "Install any missing dependencies. One of " +
          "[#{cmdline.allowed_deps.join(', ')}]",
          :type => String, :default => nil
      opt :examples, "Show some example commands with explanations"
    end
  end

  def terminal_columns
    require 'io/console'
    IO.console.winsize.last
  end

  def help_message
  <<-EOS

Transrate v#{Transrate::VERSION::STRING.dup}
by Richard Smith-Unna, Chris Boursnell, Rob Patro,
   Julian Hibberd, and Steve Kelly

DESCRIPTION:
Analyse a de-novo transcriptome assembly using three kinds of metrics:

1. sequence based (if --assembly is given)
2. read mapping based (if --left and --right are given)
3. reference based (if --reference is given)

Documentation at http://hibberdlab.com/transrate

USAGE:
transrate <options>

OPTIONS:

    EOS
  end

  def transrate_banner
    if terminal_columns > 70
      txp = '░▓▓▓^▓▓▓░'
      toptxp = txp.green
      midtxp = txp.yellow
      bottxp = txp.red
      puts <<-EOS
           _                                        _
          | |_  _ __  __ _  _ __   ___  _ __  __ _ | |_  ___
#{toptxp} | __|| '__|/ _` || '_ \\ / __|| '__|/ _` || __|/ _ \\ #{toptxp}
#{midtxp} | |_ | |  | (_| || | | |\\__ \\| |  | (_| || |_|  __/ #{midtxp}
#{bottxp}  \\__||_|   \\__,_||_| |_||___/|_|   \\__,_| \\__|\\___| #{bottxp}
      EOS
    end
    ""
  end

  def print_examples
    msg = <<-EOS

    Transrate v#{Transrate::VERSION::STRING.dup}

    EXAMPLE COMMANDS:

    # check dependencies and install any that are missing
    transrate --install-deps all

    # get the transrate score for the assembly and each contig
    transrate --assembly contigs.fa --left left.fq --right right.fq

    # basic assembly metrics only
    transrate --assembly contigs.fa

    # basic and reference-based metrics with 8 threads
    transrate --assembly contigs.fa --reference ref.fa --threads 8

    # contig and read-based metrics for two assemblies with 32 threads
    transrate --assembly one.fa,two.fa --left l.fq --right r.fq --threads 32

    EOS
    puts msg.split("\n").map{ |line| line.lstrip }.join("\n")
    exit(0)
  end

  def check_arguments
    check_dependencies
    check_loglevel
    check_assembly
    check_reference
    check_reads
    check_output
  end

  def check_loglevel
    unless %w[error info warn debug].include? @opts.loglevel
      raise TransrateError.new "Loglevel #{@opts.loglevel} is not valid. " +
      "It must be one of: error, info, warn, debug."
    end

    logger.level = Yell::Level.new @opts.loglevel.to_sym
  end

  def check_assembly
    if @opts.assembly
      @opts[:assembly] = @opts.assembly.split(',').map do |a|
        File.expand_path a
      end.join(',')
      @opts.assembly.split(',').each do |assembly_file|
        unless File.exist?(assembly_file)
          raise TransrateIOError.new "Assembly fasta file does not exist: " +
                                     " #{assembly_file}"
        end
      end
    else
      raise TransrateArgError.new "Option --assembly must be specified. " +
                                  "Try --help for help."
    end
  end

  def check_reference
    if @opts.reference
      @opts[:reference] = File.expand_path @opts.reference
      if !File.exist?(@opts.reference)
        raise TransrateIOError.new "Reference fasta file does not exist: " +
                                 " #{@opts.reference}"
      end
    end
  end

  def check_reads
    if @opts.left and @opts.right
      if @opts.left.split(",").length != @opts.right.split(",").length
        msg = "Please provide the same number of left reads as right reads"
        raise TransrateArgError.new msg
      end
      @opts[:left] = @opts.left.split(',').map { |f|
        File.expand_path f
      }.join(',')
      @opts[:right] = @opts.right.split(',').map { |f|
        File.expand_path f
      }.join(',')
      @opts.left.split(",").zip(@opts.right.split(",")).each do |left,right|
        if !File.exist?(left)
          raise TransrateIOError.new "Left read fastq file does not exist: #{left}"
        end
        if !File.exist?(right)
          raise TransrateIOError.new "Right read fastq file does not exist: #{right}"
        end
      end
    end
  end

  def check_output
    if File.exist?("#{@opts.output}/#{@outfile}")
      msg = "assemblies.csv would be overwritten in #{@opts.output}/. "
      msg << "please choose a different output directory"
      raise TransrateArgError.new msg
    end
  end

  def check_dependencies
    # Check dependencies if they are relevant to the command issued,
    # and handle any commands to install missing ones
    gem_dir = Gem.loaded_specs['transrate'].full_gem_path
    gem_deps = File.join(gem_dir, 'deps', 'deps.yaml')
    blast_dep = File.join(gem_dir, 'deps', 'blast.yaml')

    deps, read_deps, ref_deps = nil
    unless @opts.install_deps.nil?
      check_install_command

      deps = @opts.install_deps == 'all'
      read_deps = @opts.install_deps == 'read'
      ref_deps = @opts.install_deps == 'ref'
    end

    if deps || read_deps || ref_deps
      # user has requested dependency installation
      puts "Checking dependencies"
      install_missing_dependencies(deps, read_deps, ref_deps,
                                   gem_deps, blast_dep)
    else
      # no dependency installation requested, but check dependencies
      # for the commands provided are installed
      missing = []
      missing = Bindeps.missing gem_deps if @opts.left
      blast_missing = []
      blast_missing = Bindeps.missing blast_dep if @opts.reference
      print_missing_dependencies(missing, blast_missing)
    end

  end # check_dependencies

  def allowed_deps
    binkey = 'TRANSRATE_PACKAGED_BINARY'
    if ENV.has_key?(binkey) && ENV[binkey] == 'true'
      return ['ref']
    else
      return ['read', 'ref', 'all']
    end
  end

  def check_install_command
    unless allowed_deps.include? @opts.install_deps
      msg = "install-deps #{@opts.install_deps} is not valid. " +
            "You must specify one of: #{allowed_deps.join(', ')}."
      raise TransrateError.new(msg)
    end
  end

  def install_missing_dependencies(deps, read_deps, ref_deps,
                                   gem_deps, blast_dep)
    missing = []
    if deps || read_deps
      Bindeps.require gem_deps
      missing += Bindeps.missing gem_deps
    end

    if deps || ref_deps
      Bindeps.require blast_dep
      missing += Bindeps.missing blast_dep
    end

    unless missing.empty?
      list = missing.collect {|i| "#{i.name}:#{i.version}"}.join("\n - ")
      msg = "Failed to install: \n - #{list}"
      raise TransrateError.new msg
    end

    puts "All dependencies installed"
    exit
  end # install_missing_dependencies

  def print_missing_dependencies(missing, blast_missing)
    if missing.length + blast_missing.length > 0
      puts "Dependencies are missing:"

      missing.each do |dep|
        puts "  - #{dep.name} (#{dep.version})"
      end

      blast_missing.each do |dep|
        puts "  - #{dep.name} (#{dep.version})"
      end

      deps_help = {}

      deps_help[:all] =
        "To install all missing dependencies, run:\n" +
        "  transrate --install-deps all"

      deps_help[:read] =
        "To install only the read-metrics dependencies:\n" +
        "  transrate --install-deps read"

      deps_help[:read] =
        "To install only the reference-metrics dependencies:\n"
        "  transrate --install-deps ref"

      binkey = 'TRANSRATE_PACKAGED_BINARY'
      if ENV.has_key?(binkey) && ENV[binkey] == 'true'
        puts "You are running the packaged version of transrate"
        puts "This comes with the read-metrics dependencies pre-installed"
        puts ""
      end

      allowed_deps.each { |dep| puts deps_help[dep.to_sym] }
      exit 1
    end
  end

  def pretty_print_hash(hash, width, round=2)
    hash.map do |k, v|
      # show as float if there are any decimal places
      if v.to_f.round(round).to_s.split('.').last.to_i > 0
        v = v.to_f.round(round)
      end
      if v.is_a? Float
        v = v.round(round)
      end
      pad = (width - (k.to_s.length + v.to_s.length))
      pad = [pad, 0].max
      logger.info "#{k.to_s.split('_').join(' ')}" +
      "#{" " * pad}" +
      "#{v}"
    end
  end

  def concatenate_assemblies assemblies
    merged_file = @opts.merge_assemblies
    merged = {}
    assemblies.each do |file|
      Bio::FastaFormat.open(file).each do |entry|
        contig_name = "#{File.basename(file,File.extname(file))}:"
        contig_name << "#{entry.entry_id}"
        merged[contig_name] = entry.seq
      end
    end
    logger.info "Merging assemblies into one file...'#{merged_file}'"
    File.open(merged_file, "wb") do |out|
      merged.each do |name, seq|
        out.write ">#{name}\n"
        out.write "#{seq}\n"
      end
    end
    [merged_file]
  end

  def analyse_assembly(assembly, r, result_path)
    logger.info "Loading assembly: #{assembly}"
    a = Assembly.new assembly

    logger.info "Analysing assembly: #{assembly}"
    logger.info "Results will be saved in #{File.expand_path result_path}"

    contig_results = {}
    read_results = {}
    comparative_results = {}
    score, optimal, cutoff, weighted = ["NA", "NA", "NA", "NA"]

    FileUtils.mkdir_p result_path
    Dir.chdir result_path do
      transrater = Transrater.new(a, r, threads: @opts.threads)

      contig_results = contig_metrics transrater
      read_results = read_metrics transrater
      comparative_results = comparative_metrics transrater
      if (@opts.left && @opts.right)
        score, optimal, cutoff, weighted = assembly_score(assembly, transrater)
      end

      write_contig_csv a
    end

    contig_results.merge(read_results)
                  .merge(comparative_results)
                  .merge({ :assembly => assembly })
                  .merge({ :score => score })
                  .merge({ :optimal_score => optimal })
                  .merge({ :cutoff => cutoff })
                  .merge({ :weighted => weighted })

  end # analyse_assembly

  # if the user provides the same assembly multiple times it messes
  # with various conventions we use to name output files and directories
  # so we check for duplicates
  def check_assemblies_unique assembly_paths
    if assembly_paths.uniq.length != assembly_paths.length
      # find duplicates
      dupes = assembly_paths
      dupes.uniq.each do |d|
        i = dupes.index(d)
        dupes.delete_at(i) unless i.nil?
      end
      dupes.uniq!

      logger.error "There following paths were supplied more than once to --assembly:"
      dupes.each{ |d| puts "  - #{d}"}
      logger.error "Transrate uses the assembly path to set the output path"
      logger.error "so duplicates can lead to unexpected behaviour."
      logger.error "Please provide each assembly once only."
      exit 1
    end
  end

  def assembly_result_paths assemblies
    if (assemblies.length == 1)
      return [File.basename(assemblies.first, File.extname(assemblies.first))]
    end
    paths = assemblies.map { |a| File.expand_path a }
    check_assemblies_unique paths
    common_prefix = common_directory_path paths
    paths.map! { |p| p.to_path.gsub(common_prefix, "").gsub(/^\//, "") }
    paths.map { |p| assembly_result_path p }
  end

  def assembly_result_path assembly
    path = assembly.gsub(File::SEPARATOR, '_')
    File.basename(path, File.extname(path))
  end

  def common_directory_path(dirs)
    separator = File::SEPARATOR
    dir1, dir2 = dirs.minmax.map{ |dir| dir.split(separator) }
    dir1.zip(dir2).take_while{ |dn1,dn2| dn1==dn2 }.map(&:first).join(separator)
  end

  def write_assembly_csv results
    logger.info "Writing analysis results to #{@outfile}"

    CSV.open(@outfile, 'wb') do |file|

      keys = results[0].keys
      keys.delete(:assembly)
      head = [:assembly] + keys
      file << head
      results.each do |row|
        file << head.map { |x|
          entry = row[x]
          entry.is_a?(Float) ? entry.round(5) : entry
        }
      end

    end

  end # write_assembly_csv

  def contig_metrics transrater
    logger.info "Calculating contig metrics..."
    t0 = Time.now
    contig_results = transrater.assembly_metrics.basic_stats
    contig_results.merge! transrater.assembly.contig_metrics.results

    if contig_results
      logger.info "Contig metrics:"
      logger.info "-" *  @report_width
      pretty_print_hash(contig_results, @report_width)
    end

    logger.info "Contig metrics done in #{(Time.now - t0).round} seconds"
    contig_results
  end

  def read_metrics transrater
    read_results = {}
    if (@opts.left && @opts.right)
      logger.info "Calculating read diagnostics..."
      t0 = Time.now
      read_results = transrater.read_metrics(@opts.left, @opts.right).read_stats

      if read_results
        logger.info "Read mapping metrics:"
        logger.info "-" *  @report_width
        pretty_print_hash(read_results, @report_width)
      end

      logger.info "Read metrics done in #{(Time.now - t0).round} seconds"
    else
      logger.info "No reads provided, skipping read diagnostics"
    end
    read_results
  end

  def comparative_metrics transrater
    comparative_results = {}
    if @opts.reference
      logger.info "Calculating comparative metrics..."
      t0 = Time.now
      comparative_results = transrater.comparative_metrics.comp_stats

      if comparative_results
        logger.info "Comparative metrics:"
        logger.info "-" *  @report_width
        pretty_print_hash(comparative_results, @report_width)
      end

      logger.info "Comparative metrics done in #{(Time.now - t0).round} seconds"

      logger.info "-" * @report_width
    else
      logger.info "No reference provided, skipping comparative diagnostics"
    end
    comparative_results
  end

  def assembly_score(assembly, transrater)
    score = transrater.assembly_score
    weighted = transrater.weighted_score
    prefix = File.basename(assembly)
    optimal, cutoff = transrater.assembly_optimal_score prefix
    transrater.classify_contigs cutoff
    unless score.nil?
      pretty_print_hash({:TRANSRATE_ASSEMBLY_SCORE => score},
                        @report_width, 4)
      logger.info "-" * @report_width
      pretty_print_hash({:TRANSRATE_OPTIMAL_SCORE => optimal},
                        @report_width, 4)
      pretty_print_hash({:TRANSRATE_OPTIMAL_CUTOFF => cutoff},
                        @report_width, 4)
      pretty_print_hash(transrater.good_contigs, @report_width)
    end
    [score, optimal, cutoff, weighted]
  end

  def write_contig_csv a
    # write contig metrics to file for each contig
    outfile = File.expand_path "contigs.csv"
    logger.info "Writing contig metrics for each contig to #{outfile}"
    # have option to turn off, default on
    first=true
    CSV.open(outfile, 'wb') do |csv|
      a.each do |name, contig|
        basic_metrics = {:contig_name => name}.merge(contig.basic_metrics)
        if @opts.reference
          comp_metrics = contig.comparative_metrics
          basic_metrics.merge!(comp_metrics)
        end
        if @opts.left and @opts.right
          read_metrics = contig.read_metrics
          basic_metrics.merge!(read_metrics)
        end
        if first
          csv << basic_metrics.keys
          first = false
        end
        csv << basic_metrics.values.map{ |x| x.is_a?(Float) ? x.round(6) : x }
      end
    end
  end

end # Cmdline

end # Transrate
