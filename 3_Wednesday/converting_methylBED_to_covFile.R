## -------------------------------------------------------------------------
## Convert the methylBED file into something compatible with downstream analyses
## -------------------------------------------------------------------------

# All of our downstream pipelines use a .cov file from short-read data
# These .cov files are made by Bismark
# We can convert the Nanopore .methylBED file (made from modkit pileup)
# into the same format as the .cov files so it will work the same in later steps

library(readr)
library(tidyr)

# Read in data
file.list <- list.files(path = "./Nanopore_bed_files/", pattern = "*bed")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = F, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
sample <- list("F4_nanopore","F5_nanopore","F6_nanopore",
               "M1_nanopore","M2_nanopore","M3_nanopore")
names(samples) <- sample

# NOTE: after comparing WGBS and Nanopore outputs, they're all somehow
# one base out of sync... not idea how/why, probably a python/bash thing
# Add one to start and end of nanopore co-ordinates to compensate

# Edit format to make like a merged.cov file from Bismark
for(i in seq_along(samples)){
  samples[[i]] <- separate_wider_delim(samples[[i]], cols = X10, delim = " ", names_sep="unique")
  samples[[i]] <- samples[[i]][,c(1,2,3,11,12,10)]
  samples[[i]]$X2 <- samples[[i]]$X2+1
  samples[[i]]$X3 <- samples[[i]]$X3+1
  final_file <- samples[[i]]
  myfile <- file.path("./cov_files/", paste0(names(samples[i]),"_","for_methylkit.txt"))
  write.table(final_file, file=myfile, quote=F, sep="\t", row.names=F, col.names=F)
}
