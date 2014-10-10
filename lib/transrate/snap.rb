module Transrate

  class SnapError < StandardError
  end

  class Snap

    attr_reader :index_name, :sam, :read_count

    def initialize
      which_snap = Cmd.new('which snap')
      which_snap.run
      if !which_snap.status.success?
        raise SnapError.new("could not find snap in the path")
      end
      @snap = which_snap.stdout.split("\n").first

      @index_built = false
      @index_name = ""
    end

    def map_reads(file, left, right, insertsize: 200,
                  insertsd: 50, outputname: nil, threads: 8)
      raise SnapError.new("Index not built") if !@index_built
      lbase = File.basename(left.split(",").first)
      rbase = File.basename(right.split(",").first)
      index = File.basename(@index_name)
      @bam = File.expand_path("#{lbase}.#{rbase}.#{index}.bam")
      @read_count_file = "#{lbase}-#{rbase}-read_count.txt"
      realistic_dist = insertsize + (3 * insertsd)
      unless File.exists? @bam
        # construct snap command
        snapcmd = "#{@snap} paired #{@index_name}"
        left.split(",").zip(right.split(",")).each do |left, right|
          snapcmd << " #{left} #{right}"
        end
        snapcmd << " -o #{@bam}"
        snapcmd << " -I"   # ignore read IDs
        snapcmd << " -d 30" # max edit distance (function of read length?)
        snapcmd << " -t #{threads}"
        snapcmd << " -so" # sort bam output
        snapcmd << " -M"  # format cigar string
        snapcmd << " -sa" # keep all alignments, don't discard 0x100
        # run snap
        runner = Cmd.new snapcmd
        runner.run

        runner.stdout.split("\n").each do |line|
          cols = line.split(/\s+/)
          if cols.length == 13 and cols[0]=="16000"
            @read_count = cols[9].to_i / 2
            File.open("#{@read_count_file}", "wb") do |out|
              out.write("#{@read_count}\n")
            end
          end
        end

        if !runner.status.success?
          raise SnapError.new("Snap failed\n#{runner.stderr}")
        end
      else
        @read_count = 0
        if File.exist?("#{@read_count_file}")
          @read_count = File.open("#{@read_count_file}").readlines.join.to_i
        else
          left.split(",").each do |l|
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
      @bam
    end

    def build_index file, threads
      @index_name = File.basename(file, File.extname(file))
      unless File.exists?("#{@index_name}/GenomeIndexHash")
        overflow=0
        err = "Overflowed overflow table"
        while err =~ /Overflowed overflow table/ and overflow <= 950
          overflow += 50
          cmd = "#{@snap} index #{file} #{@index_name} -O#{overflow} "
          cmd << "-t#{threads}"
          runner = Cmd.new cmd
          runner.run
          err = runner.stderr
        end
        if !runner.status.success?
          msg = "Failed to build Snap index\n#{runner.stderr}"
          raise SnapError.new(msg)
        end
      end
      @index_built = true
    end

  end # Snap

end # Transrate
