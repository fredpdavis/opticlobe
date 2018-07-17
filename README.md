# opticlobe

This package contains code to analyze RNA-seq data and generate figures and
tables presented in:

A genetic, genomic, and computational resource for exploring neural circuit function
Davis FP*, Nern A*, Picard S, Reiser MB, Rubin GM, Eddy SR, Henry GL.
submitted, 2018.

Contact fred.davis@nih.gov with any questions.

## Package contents

### src

* Code, organized by language, to turn RNA-seq read files into figures.

### metadata

* text tables describing the RNA-seq samples

## Analysis overview

This section describes steps to turn RNA-seq reads to figures and tables.

## Requirements

### Software Requirements

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

Table: External software used for data processing

### Data sources

The code uses the following internal data:

1. sample descriptions. included in metadata directory
2. raw FASTQ files. available from GEO (accession GSE116969).
3. fly_neurotransmission_manuallist.txt - list of neurotransmitter-associated genes


The code also uses data from the following external data sources:

1. FlyBase: gene groups, InterPro annotation, gene ontology annotation.
2. ENSEMBL genome and transcript sequence.
3. BioMart gene ontology annotation.
4. ERCC synthetic spike-in sequence. cms_095046.txt
5. paper: Rivera-Alaba et al., 2011. Table S2.
6. paper: Takemura et al., 2013. nature12450-s3.xls
7. paper: Ozkan et al., 2013. Tables S1 and S2
8. paper: Konstantinides et al., 2018. Table S1 (cluster markers), manually entered Figure 3A cluster labels.

### Genomic data versions

1. genome assembly: BDGP6 (dm6)
2. gene models: ENSEMBL91, based on FlyBase 2017_04
3. gene groups: FlyBase 2018_02

## Processing Steps

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
