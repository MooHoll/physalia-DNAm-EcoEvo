# Bonus Goodies

Here you will find various folders with code that might be helpful for your projects. Remember, these are templates and you will need to change the parameters to suit your specific data and project goals.

### Pipelines
* RNA-Seq processing (including DNMT expression)
* SNP calling

### Extra stuff for DNA methylation
* Alternate reference genomes (incorporating SNPs)
This script takes a reference genome and substitutes SNPs of an individual into it. This helps with accurate methylation calling.
* Blast to identify DNMTs and TET
This script will show you how to identify DNMTs/TET if you're working with a non-model organism.

* [Eamonn Mallon](https://le.ac.uk/people/eamonn-mallon) has kindly agreed to share his epigenetic clock scripts [nasonia_epigenetic_ageing](https://github.com/EamonnMallon/nasonia_epigenetic_ageing). This pipeline is still a work-in-progress but an excellent starting point as it takes `differentially methylated positions` (Day 4) and uses an `elastic net regression` to identify CpGs correlated with age. There is also a cool script looking at DNA methylation entropy in the repository that you may also find useful.

* Does my organism even have DNA methylation? If you're working with a new species, you can identify from just the genome `.fa` and  annotation `.gff` whether it is likely to have a DNA methylation system or not. [Scripts here](https://github.com/MooHoll/cpg_observed_expected).
