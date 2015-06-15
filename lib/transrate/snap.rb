module Transrate

  class SnapError < TransrateError
  end

  class Snap

    require 'fix-trinity-output'
    require 'bio'

    attr_reader :index_name, :bam, :read_count

    def initialize
      which_snap = Cmd.new('which snap-aligner')
      which_snap.run
      if !which_snap.status.success?
        raise SnapError.new("could not find snap in the path")
      end
      @snap = which_snap.stdout.split("\n").first

      @index_built = false
      @index_name = ""
    end

    def build_paired_cmd l, r, threads
      cmd = "#{@snap} paired #{@index_name}"
      l.split(",").zip(r.split(",")).each do |left, right|
        cmd << " #{left} #{right}"
      end
      cmd << " -o #{@bam}"
      cmd << " -s 0 1000" # min and max distance between paired-read starts
      cmd << " -H 300000" # max seed hits to consider in paired mode
      cmd << " -h 2000" # max seed hits to consider when reverting to single
      cmd << " -d 30" # max edit distance (function of read length?)
      cmd << " -t #{threads}"
      cmd << " -b" # bind threads to cores
      cmd << " -M"  # format cigar string
      cmd << " -D 5" # extra edit distance to search. needed for -om
      cmd << " -om 5" # Output multiple alignments. extra edit distance
      cmd << " -omax 10" # max alignments per pair/read
      cmd
    end

    def map_reads(file, left, right, outputname: nil, threads: 8)
      raise SnapError.new("Index not built") if !@index_built

      lbase = File.basename(left.split(",").first)
      rbase = File.basename(right.split(",").first)
      index = File.basename(@index_name)
      @bam = File.expand_path("#{lbase}.#{rbase}.#{index}.bam")
      @read_count_file = "#{lbase}-#{rbase}-read_count.txt"

      @fixer = Fixer.new # from the fix-trinity-output gem
      unless File.exists? @bam
        snapcmd = build_paired_cmd(left, right, threads)
        runner = Cmd.new snapcmd
        runner.run
        save_readcount runner.stdout
        unless runner.status.success?
          if runner.stderr=~/Unmatched\sread\sIDs/
            logger.warn runner.stderr
            logger.warn "Unmatched read IDs. Fixing input files..."
            remap_reads(left, right, threads)
          else
            raise SnapError.new("Snap failed\n#{runner.stderr}")
          end
        end
      else
        load_readcount left
      end
      @bam
    end

    def remap_reads(left, right, threads)
      fixedleft = []
      fixedright = []
      i = 0
      left.split(",").zip(right.split(",")).each do |l, r|
        prefix = "reads-#{i}"
        @fixer.run(l, r, "#{prefix}")
        fixedleft << "#{prefix}-fixed.1.fastq"
        fixedright << "#{prefix}-fixed.2.fastq"
        i+=1
      end
      left = fixedleft.join(",")
      right = fixedright.join(",")
      File.delete(@bam)
      logger.info "Fixed input files"
      snapcmd = build_paired_cmd(left, right, threads)
      runner = Cmd.new snapcmd
      runner.run
      save_readcount runner.stdout
      unless runner.status.success?
        raise SnapError.new("Snap failed\n#{runner.stderr}")
      end
    end

    def save_readcount stdout
      stdout.split("\n").each do |line|
        cols = line.split(/\s+/)
        if cols.size > 5 and cols[0]=~/[0-9\,]+/
          @read_count = cols[0].gsub(",", "").to_i / 2
          File.open("#{@read_count_file}", "wb") do |out|
            out.write("#{@read_count}\n")
          end
        end
      end
    end

    def load_readcount reads
      @read_count = 0
      if File.exist?("#{@read_count_file}")
        @read_count = File.open("#{@read_count_file}").readlines.join.to_i
      else
        reads.split(",").each do |l|
          cmd = "wc -l #{l}"
          count = Cmd.new(cmd)
          count.run
          if count.status.success?
            @read_count += count.stdout.strip.split(/\s+/).first.to_i/4
            File.open("#{@read_count_file}", "wb") do |out|
              out.write("#{@read_count}\n")
            end
          else
            logger.warn "couldn't get number of reads from #{l}"
          end
        end
      end
    end

    def build_index file, threads
      @index_name = File.basename(file, File.extname(file))
      unless Dir.exists?(@index_name)
        cmd = "#{@snap} index #{file} #{@index_name}"
        cmd << " -s 23"
        cmd << " -t#{threads}"
        cmd << " -bSpace" # contig name terminates with space char
        runner = Cmd.new cmd
        runner.run
        if !runner.status.success?
          err = runner.stderr
          msg = "Failed to build Snap index\n#{runner.stderr}"
          raise SnapError.new(msg)
        end
      end
      @index_built = true
    end

  end # Snap

end # Transrate
