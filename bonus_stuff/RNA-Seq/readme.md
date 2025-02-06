## Differential Gene Expression from RNA-Seq

Here is a handy pipeline for differential gene expression from short-read RNA-Seq. 

There are two options for alignment and read calling, RSEM is the simplest but if you have a difficult non-model organism genome, you may prefer to use STAR directly as you can change the parameters more easily.

There is also a small script for looking at DNMT expression, assuming you know the gene IDs for your DNMTs. If not, take a look at the Blast for DNMTS scripts.