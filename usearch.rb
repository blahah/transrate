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
    subcmd = " -ublast #{query}"
    subcmd += " -db #{target}"
    subcmd += " -evalue #{evalue}"
    blast6outfile = File.basename(query) + "_" + File.basename(target) + ".b6"
    subcmd += " -userout #{blast6outfile}"
    subcmf += self.custom_output_fields
    subcmd += " -strand both"
    subcmd += " -threads #{@threads}"
    self.run subcmd
    blast6outfile
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
    self.run subcmd
  end

  def run subcmd
    subcmd += " -quiet"
    `#{@cmd}#{subcmd}`
  end

end # Usearch