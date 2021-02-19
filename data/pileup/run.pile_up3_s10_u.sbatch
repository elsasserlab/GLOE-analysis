#!/bin/bash -l

#SBATCH -A snic2020-15-9
#SBATCH -p devel
#SBATCH -n 16
#SBATCH -t 1:00:00
#SBATCH -M snowy

module load bioinfo-tools
module load samtools

module load MACS
module load deepTools
module load BEDTools

if [ $# -eq 0 ]; then
  peaks_dir="genomecov_"$(date +"%Y%m%d%H%M")
else
  peaks_dir=$1
fi

mkdir $peaks_dir
cp merge_bam_cov.sbatch "./$peaks_dir/run.$peaks_dir.sbatch"

rm -f "ngsplot.config.txt"

for filebase in GL-WT-UT-1 GL-WT-20-2 GL-WT-24hr-3 GL-KO-UT-4 GL-KO-20-5 GL-KO-24hr-6
do
#  samtools merge $filebase".hg19.fltd.bam" $filebase"-*.fltd.bam"
#  samtools index $filebase".hg19.fltd.bam"

#  bamCoverage -p 16 -e --binSize 50 --bam $filebase".hg19.fltd.bam" --normalizeUsing RPGC --effectiveGenomeSize 2864785220 -o $filebase".hg19.bs50.bw"
#  bamCoverage -p 16 --samFlagInclude 80 --bam $filebase".hg19.fltd.bam" --binSize 5 --normalizeUsing RPGC --effectiveGenomeSize 2864785220 -o $filebase".hg19.R1minus.bw"
#  bamCoverage -p 16 --samFlagInclude 64 --samFlagExclude 16 --bam $filebase".hg19.fltd.bam" --binSize 5 --normalizeUsing RPGC --effectiveGenomeSize 2864785220 -o $filebase".hg19.R1plus.bw"


#  samtools view -b -f 144 $filebase".hg19.fltd.bam" > $filebase".hg19.plus1.bam"
#  samtools index $filebase".hg19.plus1.bam"

#  samtools view -b -f 96 $filebase".hg19.fltd.bam"> $filebase".hg19.plus2.bam"
#  samtools index $filebase".hg19.plus2.bam"

#  samtools view -b -f 128 -F 16 $filebase".hg19.fltd.bam" > $filebase".hg19.minus1.bam"
#  samtools index  $filebase".hg19.minus1.bam"

#  samtools view -b -f 64 -F 32 $filebase".hg19.fltd.bam" > $filebase".hg19.minus2.bam"
#  samtools index $filebase".hg19.minus2.bam"

#  samtools merge -f $filebase".hg19.plus3.bam" $filebase".hg19.plus1.bam" $filebase".hg19.plus2.bam"
#  samtools index $filebase".hg19.plus3.bam"
#  bedtools intersect -v -abam $filebase".hg19.plus3.bam" -b consensusBlacklist.hg19.bed > $filebase".hg19.plus.bam"
#  samtools index $filebase".hg19.plus.bam"
#  samtools merge -f $filebase".hg19.minus3.bam" $filebase".hg19.minus1.bam" $filebase".hg19.minus2.bam"
#  samtools index $filebase".hg19.minus3.bam"
#  bedtools intersect -v -abam $filebase".hg19.minus3.bam" -b consensusBlacklist.hg19.bed > $filebase".hg19.minus.bam"
#  samtools index $filebase".hg19.minus.bam"

  rm -f $filebase".hg19.minus1.bam"
  rm -f $filebase".hg19.minus2.bam"
  rm -f $filebase".hg19.minus3.bam"
  rm -f $filebase".hg19.plus1.bam"
  rm -f $filebase".hg19.plus2.bam"
  rm -f $filebase".hg19.plus3.bam"
  rm -f $filebase".hg19.minus1.bam.bai"
  rm -f $filebase".hg19.minus2.bam.bai"
  rm -f $filebase".hg19.minus3.bam.bai"
  rm -f $filebase".hg19.plus1.bam.bai"
  rm -f $filebase".hg19.plus2.bam.bai"
  rm -f $filebase".hg19.plus3.bam.bai"

#  samtools view -bh -f 64 $filebase".hg19.fltd.bam" | samtools view -bh -L consensusBlacklist.hg19.bed -U $filebase".hg19.R1.bam" - > /dev/null
#  samtools index $filebase".hg19.R1.bam
  bedtools genomecov -ibam $filebase".hg19.R1.bam" -bg -strand + | awk '$4 > 2' | bedtools merge | awk '($3-$2) >= 20' | grep -v "_gl" > "./$peaks_dir/"$filebase".plus.pile.bed"
  bedtools genomecov -ibam $filebase".hg19.R1.bam" -bg -strand - | awk '$4 > 2' | bedtools merge | awk '($3-$2) >= 20' | grep -v "_gl" > "./$peaks_dir/"$filebase".minus.pile.bed"

  #bedtools merge -i "./$peaks_dir/"$filebase".plus.pile.bedgraph" > "./$peaks_dir/"$filebase".plus.merge.bed"
  #bedtools merge -i "./$peaks_dir/"$filebase".minus.pile.bedgraph" > "./$peaks_dir/"$filebase".minus.merge.bed"

  bedtools shift -i "./$peaks_dir/"$filebase".plus.pile.bed" -g hg19.chrom.sizes -s -10 > "./$peaks_dir/"$filebase".plus.pile.shifted.bed"
  bedtools shift -i "./$peaks_dir/"$filebase".minus.pile.bed" -g hg19.chrom.sizes -s 10 > "./$peaks_dir/"$filebase".minus.pile.shifted.bed"

  bedtools intersect -a "$peaks_dir/"$filebase".plus.pile.shifted.bed" -b "$peaks_dir/"$filebase".minus.pile.shifted.bed" > "./$peaks_dir/"$filebase".$peaks_dir.bed"
done


for filebase in GL-WT-UT-1 GL-WT-20-2 GL-WT-24hr-3 GL-KO-UT-4 GL-KO-20-5 GL-KO-24hr-6
do

  cp "./$peaks_dir/"$filebase".$peaks_dir.bed" "./$peaks_dir/"$filebase".$peaks_dir.uniq.bed"
  cp "./$peaks_dir/"$filebase".plus.pile.bed" "./$peaks_dir/"$filebase".plus.uniq.bed"
  cp "./$peaks_dir/"$filebase".minus.pile.bed" "./$peaks_dir/"$filebase".minus.uniq.bed"

  for sub in GL-WT-UT-1 GL-WT-20-2 GL-WT-24hr-3 GL-KO-UT-4 GL-KO-20-5 GL-KO-24hr-6
  do
    if [ $sub != $filebase ]; then
      bedtools subtract -A -a "./$peaks_dir/"$filebase".$peaks_dir.uniq.bed" -b "./$peaks_dir/"$sub".$peaks_dir.bed" > "./$peaks_dir/"$sub".$peaks_dir.temp.bed"
      cat "./$peaks_dir/"$sub".$peaks_dir.temp.bed" > "./$peaks_dir/"$filebase".$peaks_dir.uniq.bed"
      rm -f "./$peaks_dir/"$sub".$peaks_dir.temp.bed"

      bedtools subtract -A -a "./$peaks_dir/"$filebase".plus.uniq.bed" -b "./$peaks_dir/"$sub".plus.pile.bed" > "./$peaks_dir/"$filebase".plus.temp.bed"
      cat "./$peaks_dir/"$filebase".plus.temp.bed" > "./$peaks_dir/"$filebase".plus.uniq.bed"
      rm -f "./$peaks_dir/"$filebase".plus.temp.bed"

      bedtools subtract -A -a "./$peaks_dir/"$filebase".minus.uniq.bed" -b "./$peaks_dir/"$sub".minus.pile.bed" > "./$peaks_dir/"$filebase".minus.temp.bed"
      cat "./$peaks_dir/"$filebase".minus.temp.bed" > "./$peaks_dir/"$filebase".minus.uniq.bed"
      rm -f "./$peaks_dir/"$filebase".minus.temp.bed"
    fi
  done

  plus_peaks=$( cat "./$peaks_dir/"$filebase".plus.pile.bed" | wc -l )
  minus_peaks=$( cat "./$peaks_dir/"$filebase".minus.pile.bed" | wc -l )
  plus_uniq=$( cat "./$peaks_dir/"$filebase".plus.uniq.bed" | wc -l )
  minus_uniq=$( cat "./$peaks_dir/"$filebase".minus.uniq.bed" | wc -l )
  co_peaks=$( cat "./$peaks_dir/"$filebase".$peaks_dir.bed" | wc -l )
  uniq_peaks=$( cat "./$peaks_dir/"$filebase".$peaks_dir.uniq.bed" | wc -l )

  echo -e $filebase "\t" $plus_peaks "\t" $plus_uniq "\t" $minus_peaks "\t" $minus_uniq "\t" $co_peaks "\t" $uniq_peaks >> "./$peaks_dir/macs2_results.$peaks_dir.txt"

  echo -e "$filebase.hg19.fltd.bam\t./$peaks_dir/$filebase.$peaks_dir.uniq.bed\t$filebase" >> "./$peaks_dir/ngsplot.config.txt"
done


module load ngsplot

ngs.plot.r -G mm9 -R bed -T "Peaks" -C "./$peaks_dir/ngsplot.config.txt" -O "./$peaks_dir/ngsplot.uniq" -P 0 -SC global -L 1000 -IN 0 -MQ 20 -RB 0.01 -FL 150 -GO max
rm -rf "./$peaks_dir/ngsplot.uniq"
rm -f "./$peaks_dir/ngsplot.uniq.zip"
rm -f "./$peaks_dir/ngsplot.config.txt"
