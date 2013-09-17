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

  def ublast query, target, evalue="1e-5", blast6out=true
    subcmd = " -ublast #{query}"
    subcmd += " -db #{target}"
    subcmd += " -evalue #{evalue}"
    if blast6out
      blast6outfile = File.basename(query) + "_" + File.basename(target)
      subcmd += " -blast6out #{blast6outfile}"
    end
    subcmd += " -strand both"
    subcmd += " -threads #{@threads}"
    self.run subcmd
  end

  def makeudb_ublast filepath, output
    subcmd = " -makeudb_ublast #{filepath}"
    subcmd += " -output #{output}"
    self.run subcmd
  end

  def findorfs filepath, output
    subcmd = " -findorfs #{filepath}"
    subcmd += " -output #{output}"
    subcmd += " -xlat"
    subcmd += " -orfstyle 7"
    subcmd += " -threads #{@threads}"
  end

  def run subcmd
    `#{@cmd}#{subcmd}`
  end

end # Usearch