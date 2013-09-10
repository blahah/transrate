Transrate
----

Deep quality analyis and comparison of transcriptome assemblies.

## Development status

This software is in pre-alpha development and is not yet ready for deployment. 

## Analysis pipeline ##

To run the analysis pipe on a set of assemblies:

1. Run compareassemblies.rb. This will perform similarity searches and alignments, and measure the reference if one is provided.
2. Run mergeassemblymetrics.R. This will parse the metrics for all analysed assemblies into a format suitable for use in the comparison webapp.
3. Use RStudio to run the webapp.

The pipeline is under active development and in the process of being rewritten with new metrics and shaped into a coherent application.