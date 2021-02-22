#!/bin/bash -l

#SBATCH -A snic2020-15-9
#SBATCH -p devel
#SBATCH -n 20
#SBATCH -t 1:00:00

module load bioinfo-tools
module load samtools

module load MACS
module load deepTools
module load BEDTools

refdir="/proj/snic2020-6-3/bowtie2-index/"

if [ $# -eq 0 ]; then
  peaks_dir="merge_bam_cov_2M_"$(date +"%Y%m%d%H%M")
else
  peaks_dir=$1
fi

mkdir $peaks_dir
mkdir "./$peaks_dir/2M"
cp merge_bam_cov_downsample.sbatch "./$peaks_dir/"

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

  samtools view -bh -f 64 $filebase".hg19.fltd.bam" | samtools view -bh -L consensusBlacklist.hg19.bed -U $filebase".hg19.R1.bam" - > /dev/null
  samtools index $filebase".hg19.R1.bam"
  fltd=$( samtools view -c $filebase".hg19.R1.bam" )

  echo "$filebase find pileups"

  bedtools genomecov -ibam $filebase".hg19.R1.bam" -bg -strand + | awk '$4 > 2' | bedtools merge | awk '($3-$2) >= 20' | grep -v "_gl" > "./$peaks_dir/"$filebase".plus.pile.all.bed"
  bedtools genomecov -ibam $filebase".hg19.R1.bam" -bg -strand - | awk '$4 > 2' | bedtools merge | awk '($3-$2) >= 20' | grep -v "_gl" > "./$peaks_dir/"$filebase".minus.pile.all.bed"

  echo "$filebase downsample to 2M"

  samtools view -bh -f 64 $filebase".hg19.fltd.bam" | samtools view -bh -L consensusBlacklist.hg19.bed -U $filebase".hg19.R1.temp" - > /dev/null
  samtools view -h $filebase".hg19.R1.temp" | head -n 2000121 | samtools view -bh - > $filebase".hg19.R1.bam"
  rm -f $filebase".hg19.R1.temp"
  samtools index $filebase".hg19.R1.bam"

  echo "$filebase find pileups in downsampled BAM"

  bedtools genomecov -ibam $filebase".hg19.R1.bam" -bg -strand + | awk '$4 > 2' | bedtools merge | awk '($3-$2) >= 20' | grep -v "_gl" > "./$peaks_dir/"$filebase".plus.pile.bed"
  bedtools genomecov -ibam $filebase".hg19.R1.bam" -bg -strand - | awk '$4 > 2' | bedtools merge | awk '($3-$2) >= 20' | grep -v "_gl" > "./$peaks_dir/"$filebase".minus.pile.bed"

  echo "$filebase find co-peaks"

  bedtools shift -i "./$peaks_dir/"$filebase".plus.pile.bed" -g hg19.chrom.sizes -s -10 > "./$peaks_dir/"$filebase".plus.pile.shifted.bed"
  bedtools shift -i "./$peaks_dir/"$filebase".minus.pile.bed" -g hg19.chrom.sizes -s 10 > "./$peaks_dir/"$filebase".minus.pile.shifted.bed"

  bedtools intersect -a "$peaks_dir/"$filebase".plus.pile.shifted.bed" -b "$peaks_dir/"$filebase".minus.pile.shifted.bed" > "./$peaks_dir/"$filebase".copeaks.bed"

  total=$( samtools view -c -f 64 $filebase".hg19.fltd.bam" )
  fltd_2M=$( samtools view -c $filebase".hg19.R1.bam" )

  echo -e $filebase"\t"$total"\t"$fltd"\t"$fltd_2M >> "./$peaks_dir/filter_stats.txt"
done


for filebase in GL-WT-UT-1 GL-WT-20-2 GL-WT-24hr-3 GL-KO-UT-4 GL-KO-20-5 GL-KO-24hr-6
do

  echo "$filebase find unique pileups"

  cp "./$peaks_dir/"$filebase".copeaks.bed" "./$peaks_dir/"$filebase".copeaks.uniq.bed"
  cp "./$peaks_dir/"$filebase".plus.pile.bed" "./$peaks_dir/"$filebase".plus.uniq.bed"
  cp "./$peaks_dir/"$filebase".minus.pile.bed" "./$peaks_dir/"$filebase".minus.uniq.bed"

  for sub in GL-WT-UT-1 GL-WT-20-2 GL-WT-24hr-3 GL-KO-UT-4 GL-KO-20-5 GL-KO-24hr-6
  do
    if [ $sub != $filebase ]; then
      bedtools subtract -A -a "./$peaks_dir/"$filebase".copeaks.uniq.bed" -b "./$peaks_dir/"$sub".copeaks.bed" > "./$peaks_dir/"$filebase".copeaks.uniq.temp"
      cat "./$peaks_dir/"$filebase".copeaks.uniq.temp" > "./$peaks_dir/"$filebase".copeaks.uniq.bed"
      rm -f "./$peaks_dir/"$filebase".copeaks.uniq.temp"

      bedtools subtract -A -a "./$peaks_dir/"$filebase".plus.uniq.bed" -b "./$peaks_dir/"$sub".plus.pile.bed" > "./$peaks_dir/"$filebase".plus.temp.bed"
      cat "./$peaks_dir/"$filebase".plus.temp.bed" > "./$peaks_dir/"$filebase".plus.uniq.bed"
      rm -f "./$peaks_dir/"$filebase".plus.temp.bed"

      bedtools subtract -A -a "./$peaks_dir/"$filebase".minus.uniq.bed" -b "./$peaks_dir/"$sub".minus.pile.bed" > "./$peaks_dir/"$filebase".minus.temp.bed"
      cat "./$peaks_dir/"$filebase".minus.temp.bed" > "./$peaks_dir/"$filebase".minus.uniq.bed"
      rm -f "./$peaks_dir/"$filebase".minus.temp.bed"
    fi
  done

  plus_peaks_all=$( cat "./$peaks_dir/"$filebase".plus.pile.all.bed" | wc -l )
  plus_peaks=$( cat "./$peaks_dir/"$filebase".plus.pile.bed" | wc -l )
  minus_peaks_all=$( cat "./$peaks_dir/"$filebase".minus.pile.all.bed" | wc -l )
  minus_peaks=$( cat "./$peaks_dir/"$filebase".minus.pile.bed" | wc -l )
  plus_uniq=$( cat "./$peaks_dir/"$filebase".plus.uniq.bed" | wc -l )
  minus_uniq=$( cat "./$peaks_dir/"$filebase".minus.uniq.bed" | wc -l )
  co_peaks=$( cat "./$peaks_dir/"$filebase".copeaks.bed" | wc -l )
  uniq_peaks=$( cat "./$peaks_dir/"$filebase".copeaks.uniq.bed" | wc -l )

  echo -e $filebase"\t"$plus_peaks_all"\t"$plus_peaks"\t"$plus_uniq"\t"$minus_peaks_all"\t"$minus_peaks"\t"$minus_uniq"\t"$co_peaks"\t"$uniq_peaks >> "./$peaks_dir/pileup_stats.txt"

done

mkdir "./$peaks_dir/motif"

for filebase in GL-WT-UT-1 GL-WT-20-2 GL-WT-24hr-3 GL-KO-UT-4 GL-KO-20-5 GL-KO-24hr-6
do

  echo "$filebase extract read sequences for motif"

  samtools view -bL "./$peaks_dir/"$filebase".plus.uniq.bed" $filebase".hg19.R1.bam" | bedtools bamtobed > "./$peaks_dir/"$filebase".plus.reads.bed"
  samtools view -bL "./$peaks_dir/"$filebase".minus.uniq.bed" $filebase".hg19.R1.bam" | bedtools bamtobed > "./$peaks_dir/"$filebase".minus.reads.bed"

  cat "./$peaks_dir/"$filebase".plus.reads.bed" "./$peaks_dir/"$filebase".minus.reads.bed" > "./$peaks_dir/"$filebase".pileup.reads.bed"
  rm -f "./$peaks_dir/"$filebase".plus.reads.bed"
  rm -f "./$peaks_dir/"$filebase".minus.reads.bed"

  bedtools shift -i "./$peaks_dir/"$filebase".pileup.reads.bed" -g hg19.chrom.sizes -p -20 -m 20 > "./$peaks_dir/motif/"$filebase".pileup.motif.range.bed"
  bedtools getfasta -fi $refdir"hg19.fa" -bed "./$peaks_dir/motif/"$filebase".pileup.motif.range.bed" -s | tr '[:lower:]' '[:upper:]' > "./$peaks_dir/motif/"$filebase".pileup.motif.reads.fasta"

  grep -v ">" "./$peaks_dir/motif/"$filebase".pileup.motif.reads.fasta" | wc -l > "./$peaks_dir/motif/"$filebase".count.txt"
  grep -v ">" "./$peaks_dir/motif/"$filebase".pileup.motif.reads.fasta" | cut -c21-26 | sort | uniq -c | sort -nr > "./$peaks_dir/motif/"$filebase".motif_shift20_counts_pos1-6.txt"
  grep -v ">" "./$peaks_dir/motif/"$filebase".pileup.motif.reads.fasta" | cut -c23-26 | sort | uniq -c | sort -nr > "./$peaks_dir/motif/"$filebase"motif_shift20_counts_pos3-6.txt"
  grep -v ">" "./$peaks_dir/motif/"$filebase".pileup.motif.reads.fasta" | cut -c24-25 | sort | uniq -c | sort -nr > "./$peaks_dir/motif/"$filebase"motif_shift20_counts_pos4-5.txt"
  grep -v ">" "./$peaks_dir/motif/"$filebase".pileup.motif.reads.fasta" | cut -c23-25 | sort | uniq -c | sort -nr > "./$peaks_dir/motif/"$filebase"motif_shift20_counts_pos3-5.txt"

done

mkdir "./$peaks_dir/pileup/"
for filebase1 in GL-WT-UT-1 GL-WT-20-2 GL-WT-24hr-3 GL-KO-UT-4 GL-KO-20-5 GL-KO-24hr-6
do

  for rep in 1 2 3
  do
    filebase=$filebase1"-"$rep

    echo "$filebase reads under pileup fraction"

#    samtools view -bh -f 64 $filebase".hg19.fltd.bam" | samtools view -U $filebase".hg19.R1.bam"  -L consensusBlacklist.hg19.bed - > /dev/null
#    samtools index $filebase".hg19.R1.bam"
#    fwd=$( samtools view -cL "./$peaks_dir/"$filebase1".plus.pile.bed" $filebase".hg19.R1.bam" )
#    rev=$( samtools view -cL "./$peaks_dir/"$filebase1".minus.pile.bed" $filebase".hg19.R1.bam" )
#    fltd=$( samtools view -c $filebase".hg19.R1.bam" )
    rm -f $filebase".hg19.R1.bam"
    rm -f $filebase".hg19.R1.bam.bai"

    fwd=$( samtools view -bh -f 64 $filebase".hg19.fltd.bam" | bedtools intersect -v -abam stdin -b consensusBlacklist.hg19.bed | samtools view -cL "./$peaks_dir/"$filebase1".plus.pile.bed" - )
    rev=$( samtools view -bh -f 64 $filebase".hg19.fltd.bam" | bedtools intersect -v -abam stdin -b consensusBlacklist.hg19.bed | samtools view -cL "./$peaks_dir/"$filebase1".minus.pile.bed" - )

    echo "$filebase reads under piles:" $fwd $rev
    sum=`expr $fwd + $rev`

    total=$( samtools view -c -f 64 $filebase".hg19.fltd.bam" )
    fltd=$( samtools view -bh -f 64 $filebase".hg19.fltd.bam" | bedtools intersect -v -abam stdin -b consensusBlacklist.hg19.bed | samtools view -c )

    echo -e $filebase"\t"$total"\t"$fltd"\t"$sum >> "./$peaks_dir/pileup/pileup_fractions.txt"

  done

done
