#!/bin/bash -l

#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 1:00:00
#SBATCH -M snowy

module load bioinfo-tools
module load samtools

for filebase in GL-WT-UT-1 GL-WT-20-2 GL-WT-24hr-3 GL-KO-UT-4 GL-KO-20-5 GL-KO-24hr-6
do
  fwd=$( samtools view -f 64 -cL "./pile_up3_s10_u_final/"$filebase".plus.pile.bed" $filebase".hg19.R1.bam" )
  rev=$( samtools view  -f 64 -cL "./pile_up3_s10_u_final/"$filebase".minus.pile.bed" $filebase".hg19.R1.bam" )
  sum=`expr $fwd + $rev`

  total=$( samtools view -c -f 64 $filebase".hg19.fltd.bam" )

  echo -e $filebase "\t" $sum "\t" $total >> pile_up3_fraction_merged.txt
done
