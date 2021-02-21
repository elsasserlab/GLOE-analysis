#!/bin/bash -l

#SBATCH -A snic2020-15-9
#SBATCH -p devel
#SBATCH -n 16
#SBATCH -t 1:00:00
#SBATCH -M snowy

module load bioinfo-tools
module load samtools
module load BEDTools

out="count_reads_under_pileup/"
mkdir $out

for filebase1 in GL-WT-UT-1 GL-WT-20-2 GL-WT-24hr-3 GL-KO-UT-4 GL-KO-20-5 GL-KO-24hr-6
do
  for rep in 1 2 3
  do
    filebase=$filebase1"-"$rep

    samtools view -bh -f 64 $filebase".hg19.fltd.bam" | bedtools intersect -v -abam stdin -b consensusBlacklist.hg19.bed | bedtools genomecov -bg -strand + -ibam stdin | awk '$4 > 2' | bedtools merge | awk '($3-$2) >= 20' | grep -v "_gl" > "./$out/"$filebase".plus.pile3.bed"
    samtools view -bh -f 64 $filebase".hg19.fltd.bam" | bedtools intersect -v -abam stdin -b consensusBlacklist.hg19.bed | bedtools genomecov -bg -strand - -ibam stdin | awk '$4 > 2' | bedtools merge | awk '($3-$2) >= 20' | grep -v "_gl" > "./$out/"$filebase".minus.pile3.bed"

    fwd=$( samtools view -f 64 -cL "./$out/"$filebase".plus.pile3.bed" $filebase".hg19.fltd.bam" )
    rev=$( samtools view -f 64 -cL "./$out/"$filebase".minus.pile3.bed" $filebase".hg19.fltd.bam" )
    sum=`expr $fwd + $rev`

    total=$( samtools view -c -f 64  $filebase".hg19.fltd.bam" )
    fltd=$( samtools view -bh -f 64 $filebase".hg19.fltd.bam" | bedtools intersect -v -abam stdin -b consensusBlacklist.hg19.bed | wc -l )

    echo -e $filebase"\t"$total"\t"$fltd"\t"$sum >> "./$out/pile_up3_fraction.txt"

    rm -f "./$out/"$filebase".plus.pile3.bed"
    rm -f "./$out/"$filebase".minus.pile3.bed"
  done
done
