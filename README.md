# opticlobe

This package contains code to analyze RNA-seq data and generate figures and
tables presented in:

A genetic, genomic, and computational resource for exploring neural circuit function  
Davis FP*, Nern A*, Picard S, Reiser MB, Rubin GM, Eddy SR, Henry GL.  
[bioRxiv 2018](http://doi.org/10.1101/385476)

Version: updated Aug 27, 2019 to reflect latest preprint version

Contact fredpdavis@gmail.com with any questions.

## Package contents

- src - code, organized by language, to turn RNA-seq read files into figures.
- metadata - text tables describing RNA-seq samples
- data - text data files used by code

## Requirements

### Data

The data used by this package comes from several sources. We include nearly all
data files expected by the R and LSF shell programs, along with README files
describing the contents. The only exceptions are large files (eg, genome
sequence, gene annotations, transcript sequences, RNA-seq alignment indices),
which we do not provide but describe in README files how to obtain or build.

| Data                    | Source                                                                                                                 |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| RNA-seq FASTQ files     | [GEO accession GSE116969](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116969)                                |
| genome sequence         | [ENSEMBL release 91; based on BDGP6 (dm6)](http://dec2017.archive.ensembl.org/Drosophila_melanogaster/Info/Index)      |
| gene structures         | [ENSEMBL release 91; based on FlyBase 2017_04](http://dec2017.archive.ensembl.org/Drosophila_melanogaster/Info/Index)  |
| FlyBase gene groups     | [FlyBase 2018_02](http://fb2018_02.flybase.org/)


We also used data from the following papers:

- [Rivera-Alaba et al., Curr Biol 2011.](http://dx.doi.org/10.1016/j.cub.2011.10.022)
- [Takemura et al., Nature 2013.](http://dx.doi.org/10.1038/nature12450)
- [Ozkan et al., Cell 2013.](http://dx.doi.org/10.1016/j.cell.2013.06.006)
- [Konstantinides et al., Cell 2018.](http://dx.doi.org/10.1016/j.cell.2018.05.021)
- [Davie et al., Cell 2018.](http://dx.doi.org/10.1016/j.cell.2018.05.057)


### Software

We used the following software on an LSF linux compute cluster.

|  #  | Software               |  Source                                    |
| --- | ---------------------- | ------------------------------------------ |
|  1  | seqtk (2012 Oct 16)    | https://github.com/lh3/seqtk               |
|  2  | kallisto 0.43.1        | https://pachterlab.github.io/kallisto/     |
|  3  | STAR 2.5.3c            | https://github.com/alexdobin/STAR          |
|  4  | bedtools 2.15.0        | http://bedtools.readthedocs.io/            |
|  5  | samtools 0.1.19        | https://github.com/samtools/samtools       |
|  6  | UCSC bedGraphToBigWig  | http://hgdownload.cse.ucsc.edu/admin/exe   |
|  7  | picard 1.9.1           | http://broadinstitute.github.io/picard/    |
|  8  | R v3.3.1               | https://www.r-project.org/                 |
|  9  | RStan/Stan             | http://mc-stan.org/users/interfaces/rstan  |

## Analysis overview

The analysis includes three parts: (1) RNA-seq read processing, (2) modeling
gene expression states, and (3) generating figures and tables.

## 1. RNA-seq read processing

This step is implemented in an LSF script that trims the reads (seqtk),
pseudo-aligns them to the transcriptome (kallisto) and aligns them to the
genome (STAR), evaluates quality metrics (picard), and generates bigwig
tracks for visualization (samtools, ucsc kent tools)

```sh
bsub < process_sample.lsf.sh
```

## 2. Modeling gene expression states

We implemented this step in an R routine that uses the RStan interface to 
the Stan statistical inference engine to model gene expression states.
This step uses an LSF compute cluster for parallelization.

The routine assume 200 cluster nodes with 4 CPUs each (values specified in the
setSpecs() routine). On our cluster, this step takes ~ 3 hours in total.

```R
source("../../src/R/analyzeOpticLobeExpr.R")
dat <- main(returnData = TRUE)
runModel(dat, mode="submit")
dat$pFit <- runModel(dat, mode="merge")
dat <- main(dat,returnData = TRUE)
saveRDS(dat, paste0(dat$specs$outDir,"/full_dataset.rds"))
```

## 3. Generating figures and tables

An R routine generates all the figures in final form, except for a few where we
manually repositioned labels for clarity 

```R
makeFigs(dat)
```
