# Transcriptome analysis tools #

Use this repo to co-ordinate work in progress (i.e. pre-publication) on transcriptome analysis tools.

## Analysis pipeline ##

To run the analysis pipe on a set of assemblies:

1. Run compareassemblies.rb. This will perform similarity searches and alignments, and measure the reference if one is provided.
2. Run mergeassemblymetrics.R. This will parse the metrics for all analysed assemblies into a format suitable for use in the comparison webapp.
3. Use RStudio to run the webapp.

The pipeline is under active development and will ultimately be provided with a single script to run the entire pipeline and optionally launch the