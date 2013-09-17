class Usearch

  DB_EXT = '.udb'

  require 'which'
  include Which

  def initialize
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
    self.run subcmd
  end

  def makeudb_ublast filepath, output
    subcmd = " -makeudb_ublast #{filepath}"
    subcmd += " -output #{output}"
    self.run subcmd
  end

  def run subcmd
    `#{@cmd}#{subcmd}`
  end

end # Usearch