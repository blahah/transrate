module Transrate

  class Usearch

    require 'which'
    include Which

    def initialize threads=8
      @threads = threads
      paths = which('usearch')
      if paths.empty?
        raise "usearch not found in path. Please ensure usearch is installed and aliased as 'usearch' in your path."
      end
      @cmd = paths.first
    end

    def custom_output_fields
      " -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+tcov"
    end

    def ublast query, target, evalue="1e-5"
      blast6outfile = "#{File.basename(query)}_#{File.basename(target)}.b6"
      unless File.exists? blast6outfile
      subcmd = " -ublast #{query}"
        subcmd += " -db #{target}"
        subcmd += " -evalue #{evalue}"
        subcmd += " -userout #{blast6outfile}"
        subcmd += self.custom_output_fields
        subcmd += " -strand both"
        subcmd += " -threads #{@threads}"
        self.run subcmd
      end
      blast6outfile
    end

    def makeudb_ublast filepath, output
      unless File.exists? output
        subcmd = " -makeudb_ublast #{filepath}"
        subcmd += " -output #{output}"
        self.run subcmd
      end
    end

    def findorfs filepath, output
      unless File.exists? output
        subcmd = " -findorfs #{filepath}"
        subcmd += " -output #{output}"
        subcmd += " -xlat"
        subcmd += " -orfstyle 7"
        self.run subcmd
      end
    end

    def run subcmd
      subcmd += " -quiet"
      `#{@cmd}#{subcmd} 2&>1`
    end

  end # Usearch

end # Transrate
