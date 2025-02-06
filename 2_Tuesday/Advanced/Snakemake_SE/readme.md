## Snakemake for workflow efficency

Snakemake is a great way to string together a bunch of commands to improve efficency of your pipelines.

[Documentation here.](https://snakemake.readthedocs.io/en/stable/).

The files provided in this directory are a good starting point for setting up work own workflow, but will depend on your cluster.

Here I use a conda environment to install all of the packages I will use, e.g. Bismark, Samtools etc. Then I can submit a job directly to my SLURM based cluster.
