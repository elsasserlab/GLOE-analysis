#!/bin/bash -l

#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 8:00:00

bedpath="/proj/snic2020-6-3/bed_mm9/"
indexpath="/proj/snic2020-6-3/bowtie2-index"

filebase=$1

module load bioinfo-tools
module load bowtie2
module load samtools
module load IGVtools
module load FastQC
module load picard
module load deepTools



#mkdir "pipeline"
#mkdir "pipeline/$2"

R1=$filebase"_R1.fastq.gz"
R2=$filebase"_R2.fastq.gz"

#if [ ! -f $filebase".hg19.all.bam" ]
#then
bowtie2 -p 8 -x $indexpath/mm9 -q -1 $R1 -2 $R2 --very-fast --un-conc-gz $filebase"_unR%.fastq.gz" 2> $filebase".mm9.all.bowtie.log" > /dev/null

bowtie2 -p 8 -x $indexpath/hg19 -q -1 $filebase"_unR1.fastq.gz" -2 $filebase"_unR2.fastq.gz" --fast 2> $filebase".hg19.all.bowtie.log" | samtools view -bS -F 4 - | samtools sort -o $filebase".hg19.all.bam" -
samtools index $filebase".hg19.all.bam"
samtools flagstat $filebase".hg19.all.bam" > $filebase".hg19.all.flagstat.txt"
samtools idxstats $filebase".hg19.all.bam" > $filebase".hg19.all.idxstats.txt"
#fi

#if [ ! -f $filebase".hg19.fltd.bam" ]
#then
java -jar $PICARD_HOME/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE I=$filebase".hg19.all.bam" O=$filebase".hg19.fltd.bam" M=$filebase".hg19.all.MarkDuplicates.txt"
samtools index $filebase".hg19.fltd.bam"
samtools flagstat $filebase".hg19.fltd.bam" > $filebase".hg19.fltd.flagstat.txt"
samtools idxstats $filebase".hg19.fltd.bam" > $filebase".hg19.fltd.idxstats.txt"
#fi

#samtools view -F 16 -b -o $filebase".hg19.plus.bam" $filebase".hg19.fltd.bam"
#samtools view -f 16 -b -o $filebase".hg19.minus.bam" $filebase".hg19.fltd.bam"

#samtools index $filebase".hg19.plus.bam"
#samtools flagstat $filebase".hg19.plus.bam" > $filebase".hg19.plus.flagstat.txt"
#samtools idxstats $filebase".hg19.plus.bam" > $filebase".hg19.plus.idxstats.txt"

#samtools index $filebase".hg19.minus.bam"
#samtools flagstat $filebase".hg19.minus.bam" > $filebase".hg19.minus.flagstat.txt"
#samtools idxstats $filebase".hg19.minus.bam" > $filebase".hg19.minus.idxstats.txt"

#if [ ! -f $filebase".hg19.bw" ]
#then
java -jar $PICARD_HOME/picard.jar CollectInsertSizeMetrics I=$filebase".hg19.fltd.bam" O=$filebase".hg19.insert_sizes.txt" H=$filebase".hg19.insert_sizes.pdf" M=0.5 &
#igvtools count --stands first $filebase".hg19.fltd.bam" $filebase".hg19.tdf" hg19.chrom.sizes &
#bamCoverage --bam $filebase".hg19.fltd.bam" --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2864785220 -o $filebase".hg19.bs1.bw" &
#bamCoverage -e --bam $filebase".hg19.fltd.bam" --binSize 5 --normalizeUsing RPGC --effectiveGenomeSize 2864785220 -o $filebase".hg19.ext.bs5.bw" &
bamCoverage -p 2 -e --binSize 50 --bam $filebase".hg19.fltd.bam" --normalizeUsing RPGC --effectiveGenomeSize 2864785220 -o $filebase".hg19.bs50.bw" &
#bamCoverage -e --minFragmentLength 30 --maxFragmentLength 100 --bam $filebase".hg19.fltd.bam" --normalizeUsing RPGC --effectiveGenomeSize 2864785220 -o $filebase".hg19.small.bw" &
#bamCoverage -e --minFragmentLength 150 --maxFragmentLength 300  --bam $filebase".hg19.fltd.bam" --normalizeUsing RPGC --effectiveGenomeSize 2864785220 -o $filebase".hg19.nucl.bw" &
#fi
bamCoverage -p 3 --samFlagInclude 80 --bam $filebase".hg19.fltd.bam" --binSize 5 --normalizeUsing RPGC --effectiveGenomeSize 2864785220 -o $filebase".hg19.R1minus.bw" &
bamCoverage -p 3 --samFlagInclude 64 --samFlagExclude 16 --bam $filebase".hg19.fltd.bam" --binSize 5 --normalizeUsing RPGC --effectiveGenomeSize 2864785220 -o $filebase".hg19.R1plus.bw" &
wait
