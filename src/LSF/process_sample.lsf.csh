#!/bin/csh
#BSUB -o lsfOutDir/%J.%I.out
#BSUB -e lsfOutDir/%J.%I.err
#BSUB -n 4
#BSUB -W 10:00
#BSUB -J process_samples[1-266]

source /misc/lsf/conf/cshrc.lsf
set NUMCPU=4

set origdir=`pwd`
set curhost=`hostname`

set scratchdir="/scratch/davisf/nodejobtmp-$$-$LSB_JOBID-$LSB_JOBINDEX"
mkdir -p $scratchdir
cd $scratchdir


## SET FILE PATHS AND OPTIONS

set BASEDIR="/groups/eddy/home/davisf/work/genomics/eddy_seq_data_analysis/opticlobe"
set PERLDIR="$BASEDIR/src/perl"
set SWDIR="/groups/eddy/home/davisf/software"
set SAMPLEINFO_FN="$BASEDIR/metadata/opticlobe_ms_rnaseq_samples.txt"

set fastqs=( `perl $PERLDIR/txt2tasklist.pl $SAMPLEINFO_FN fastq` )
set sampleids=( `perl $PERLDIR/txt2tasklist.pl $SAMPLEINFO_FN sampleID` )

set FASTQ=$fastqs[$LSB_JOBINDEX]
set SAMPLEID=$sampleids[$LSB_JOBINDEX]

set FASTQ_FN="$BASEDIR/data/fastq/$FASTQ"

## STEP OPTIONS / FILE/PATH specification

### BINARY PATHS

set SEQTK_BIN="$SWDIR/seqtk/seqtk"
set KALLISTO_BIN="${SWDIR}/kallisto/kallisto_linux-v0.43.1/kallisto"
set STAR_BIN="$SWDIR/STAR/STAR_2.5.3a/STAR"
set SAMTOOLS_BIN="$SWDIR/samtools/samtools-0.1.19/samtools"
set BEDTOOLS_DIR="$SWDIR/bedtools/BEDTools-Version-2.15.0/bin"
set WIGTOBIGWIG_BIN="$SWDIR/bigWig/wigToBigWig"
set COUNTFASTQ_BIN="$PERLDIR/count_fastq_reads.pl"
set SCALEBED_BIN="$PERLDIR/scale_bed_rpm.pl"
set CONVERT_ENSEMBL_TO_UCSC_BED_BIN="$PERLDIR/convert_ensembl_to_ucsc_chrom.pl"
set PICARDSTATS_BIN="/usr/local/java/bin/java -jar $SWDIR/picard/picard-tools-1.91/CollectRnaSeqMetrics.jar"


### INPUT FILES
set KALLISTO_INDEX="${BASEDIR}/data/kallisto_files.BDGP6.91/BDGP6.91.ERCC.INTACT.ncrna.kallisto_index"
set STAR_FILE_DIR="$BASEDIR/data/star_files.BDGP6.91"
set STAR_INDEX="$STAR_FILE_DIR/BDGP6.91.ERCC.INTACT.star_index"
set RRNA_BED_FN="$STAR_FILE_DIR/rRNA_regular_chromosomes.bed"
set NONRRNA_BED_FN="$STAR_FILE_DIR/regular_chromosomes.nonribosomal_regions.bed"
set REFFLAT_FN="$BASEDIR/data/picard_files.BDGP6.91/refFlat.noERCC.txt.gz"
set ERCC_REFFLAT_FN="$BASEDIR/data/picard_files.BDGP6.91/refFlat.ERCConly.txt.gz"
set GENOME_FASTAF="$STAR_FILE_DIR/BDGP6.91.ERCC.INTACT.fa"


### OUTPUT DIRECTORIES
set TRIM_OUTDIR="$BASEDIR/results/RNAseq/trim_reads/$SAMPLEID"
set KALLISTO_OUTDIR="${BASEDIR}/results/RNAseq/kallisto.BDGP6.91/${SAMPLEID}"
set STAR_OUTDIR="$BASEDIR/results/RNAseq/star_align.BDGP6.91/$SAMPLEID"
set PICARDSTATS_OUTDIR="$BASEDIR/results/RNAseq/picard_stats.star.BDGP6.91/$SAMPLEID"

### OUTPUT FILES

set TRIM_FASTQF="$TRIM_OUTDIR/${SAMPLEID}_trim.fastq.gz"

set STAR_OUT_SAM_FN="${STAR_OUTDIR}/Aligned.out.sam"
set STAR_OUT_BAM_PREFIX="${STAR_OUTDIR}/${SAMPLEID}.sorted"
set STAR_OUT_BAM_FN="${STAR_OUT_BAM_PREFIX}.bam"
set STAR_OUT_BED_FN="${STAR_OUTDIR}/${SAMPLEID}.star_scaled10M_nonRRNA.bed"
set STAR_OUT_BW_FN="$STAR_OUTDIR/${SAMPLEID}.star_scaled10M_nonRRNA.bw"

set TEMP_CHROMSIZES_FN="$STAR_OUTDIR/bam_chromsizes.txt"
set TEMP_CHROMSIZES_UCSC_FN="$STAR_OUTDIR/bam_chromsizes_ucsc.txt"
set ALICOUNT_FN="$STAR_OUTDIR/$SAMPLEID.alignment_stats.txt"

set PICARDSTATS_TXT_FN="${PICARDSTATS_OUTDIR}/${SAMPLEID}.noERCC.picard_rnaseq_report.txt"
set PICARDSTATS_PDF_FN="${PICARDSTATS_OUTDIR}/${SAMPLEID}.noERCC.picard_rnaseq_report.pdf"
set PICARDSTATS_ERCC_TXT_FN="${PICARDSTATS_OUTDIR}/${SAMPLEID}.ERCConly.picard_rnaseq_report.txt"
set PICARDSTATS_ERCC_PDF_FN="${PICARDSTATS_OUTDIR}/${SAMPLEID}.ERCConly.picard_rnaseq_report.pdf"

## STEP 0. MAKE SURE FASTQ EXISTS

if (! -e $FASTQ_FN ) then
   echo "FATAL ERROR: FASTQ file not found: $FASTQ_FN"
   exit;
endif


### STEP 1. TRIM
set SEQTK_OPTIONS="trimfq -b 5 ${FASTQ_FN}"

### STEP 2. KALLISTO PSEUDOALIGN

set KALLISTO_OPTIONS="quant -i $KALLISTO_INDEX --single -l 250 -s 50 -b 50 -t $NUMCPU -o $KALLISTO_OUTDIR $TRIM_FASTQF"

### STEP 3A. GENOME ALIGN, QC, TRACK VIZ

set STAR_OPTIONS="--genomeDir $STAR_INDEX --readFilesIn $TRIM_FASTQF --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterType BySJout --runThreadN $NUMCPU --outFileNamePrefix $STAR_OUTDIR/ --readFilesCommand zcat"


set PICARDSTATS_OPTIONS1="REF_FLAT=$REFFLAT_FN STRAND_SPECIFICITY=NONE REFERENCE_SEQUENCE=$GENOME_FASTAF INPUT=$STAR_OUT_BAM_FN CHART_OUTPUT=$PICARDSTATS_PDF_FN OUTPUT=$PICARDSTATS_TXT_FN"

set PICARDSTATS_OPTIONS2="REF_FLAT=$ERCC_REFFLAT_FN STRAND_SPECIFICITY=NONE REFERENCE_SEQUENCE=$GENOME_FASTAF INPUT=$STAR_OUT_BAM_FN CHART_OUTPUT=$PICARDSTATS_ERCC_PDF_FN OUTPUT=$PICARDSTATS_ERCC_TXT_FN"

## SETUP NECESSARY DIRECTORIES

foreach T_OUTDIR ( $TRIM_OUTDIR $KALLISTO_OUTDIR $STAR_OUTDIR $PICARDSTATS_OUTDIR )
   if (! -e $T_OUTDIR) then
      mkdir -p $T_OUTDIR
   endif
end



set curtime=`date`
echo "# cluster run started on $curhost at $curtime"
echo "# Processing: $SAMPLEID ($FASTQ)" ;

set echo

## STEP 1. READ TRIMMING
set curtime=`date`
echo "# STEP 1. READ TRIMMING ($curtime)"

$SEQTK_BIN $SEQTK_OPTIONS | gzip > ${TRIM_FASTQF}


## STEP 2. KALLISTO

set curtime=`date`
echo "# STEP 2. kallisto ($curtime)"

$KALLISTO_BIN $KALLISTO_OPTIONS


## STEP 3. STAR

set curtime=`date`
echo "# STEP 3. STAR ($curtime)"

$STAR_BIN $STAR_OPTIONS
$SAMTOOLS_BIN view -@ $NUMCPU -b -S $STAR_OUT_SAM_FN | $SAMTOOLS_BIN sort -@ $NUMCPU -m 5G - $STAR_OUT_BAM_PREFIX


## STEP 3b. alignment stats

set curtime=`date`
echo "# STEP 3b. Counting alignment statistics ($curtime)"

set NUMREADS=`perl $COUNTFASTQ_BIN $TRIM_FASTQF`
echo "# Number of reads: $NUMREADS" > $ALICOUNT_FN

set NUMALN_READS=`$SAMTOOLS_BIN view $STAR_OUT_BAM_FN | cut -f1 | sort -u -T${scratchdir} | wc -l`
echo "# Number of aligned reads: $NUMALN_READS" >> $ALICOUNT_FN

set NUMALN_RRNA_READS=`$BEDTOOLS_DIR/intersectBed -abam $STAR_OUT_BAM_FN -b $RRNA_BED_FN | $SAMTOOLS_BIN view - | cut -f1 | sort -u -T${scratchdir} | wc -l`
echo "# Number of reads aligned to rRNA loci: $NUMALN_RRNA_READS" >> $ALICOUNT_FN

set NUMALNS=`$SAMTOOLS_BIN flagstat $STAR_OUT_BAM_FN | grep -m 1 mapped | awk '{print $1}'`
echo "# Number of alignments: $NUMALNS" >> $ALICOUNT_FN

set NUMALNS_NONRRNA=`$BEDTOOLS_DIR/intersectBed -abam $STAR_OUT_BAM_FN -b $NONRRNA_BED_FN | $SAMTOOLS_BIN view -c -`
echo "# Number of non-ribosomal, non-ERCC, non-INTACT alignments: $NUMALNS_NONRRNA" >> $ALICOUNT_FN

## STEP 3c. make bigwig tracks

set curtime=`date`
echo "# STEP 3c. Making bigwig tracks ($curtime)"

$SAMTOOLS_BIN view -H $STAR_OUT_BAM_FN | grep "@SQ" | sed -e "s/SN://" -e "s/LN://" | cut -f 2,3 > $TEMP_CHROMSIZES_FN
perl $CONVERT_ENSEMBL_TO_UCSC_BED_BIN < $TEMP_CHROMSIZES_FN > $TEMP_CHROMSIZES_UCSC_FN

${BEDTOOLS_DIR}/genomeCoverageBed -split -bg -ibam $STAR_OUT_BAM_FN -g ${TEMP_CHROMSIZES_FN} | perl $SCALEBED_BIN $NUMALNS_NONRRNA 10000000 > $STAR_OUT_BED_FN
cat $STAR_OUT_BED_FN | perl $CONVERT_ENSEMBL_TO_UCSC_BED_BIN | $WIGTOBIGWIG_BIN -clip stdin $TEMP_CHROMSIZES_UCSC_FN $STAR_OUT_BW_FN 


## STEP 4. PICARD QC

set curtime=`date`
echo "# STEP 4. Running picard ($curtime)"

$PICARDSTATS_BIN $PICARDSTATS_OPTIONS1
$PICARDSTATS_BIN $PICARDSTATS_OPTIONS2


## STEP 5. CLEANUP

set curtime=`date`
echo "# STEP 5. Cleaning up ($curtime)"

cd $origdir
rm $TRIM_FASTQF
rmdir $TRIM_OUTDIR
rm $STAR_OUT_SAM_FN
rm $STAR_OUT_BAM_FN
rm $STAR_OUT_BED_FN
rm -rf $scratchdir

set curtime=`date`
echo "#cluster run finished on $curhost at $curtime"
