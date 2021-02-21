
pileup_regions=$1
R1_bam=$2
ref="/Volumes/DATA/DATA/genome/"



samtools view -bL $pileup_regions".plus.uniq.bed" $R1_bam | bamtobed > $pileup_regions".plus.reads.bed"
samtools view -bL $pileup_regions".minus.uniq.bed" $R1_bam | bamtobed > $pileup_regions".minus.reads.bed"
cat $pileup_regions".plus.reads.bed" $pileup_regions".minus.reads.bed" > $pileup_regions".pileup.reads.bed"

bedtools shift -i $pileup_regions".pileup.reads.bed" -g $ref"hg19.chrom.sizes" -p -20 -m 20 > $pileup_regions".pileup.reads.shifted.bed"

bedtools getfasta -fi $ref"hg19.fa" -bed $pileup_regions".pileup.reads.shifted.bed" -s | tr '[:lower:]' '[:upper:]' > $pileup_regions".pileup.reads.shifted.fasta"

grep -v ">" $pileup_regions".pileup.reads.shifted.fasta" | wc -l > $pileup_regions".count.txt"
grep -v ">" $pileup_regions".pileup.reads.shifted.fasta" | cut -c21-26 | sort | uniq -c | sort > $pileup_regions".motif_shift20_counts_pos1-6.txt"
grep -v ">" $pileup_regions".pileup.reads.shifted.fasta" | cut -c23-26 | sort | uniq -c | sort > $pileup_regions"motif_shift20_counts_pos3-6.txt"
grep -v ">" $pileup_regions".pileup.reads.shifted.fasta" | cut -c24-25 | sort | uniq -c | sort > $pileup_regions"motif_shift20_counts_pos4-5.txt"
grep -v ">" $pileup_regions".pileup.reads.shifted.fasta" | cut -c23-25 | sort | uniq -c | sort > $pileup_regions"motif_shift20_counts_pos3-5.txt"

# not implemented is padding of fasta s/^([catgn]{50}).+\n/$1\n/

weblogo -F pdf -n 50 < $pileup_regions".pileup.reads.shifted.fasta" > $pileup_regions".pileup.reads.logo.pdf"
