#!/bin/bash


	basepath="."


	ls "$basepath/fastq/"*_R1.fastq.gz | sed -e "s/_R1.fastq.gz//g" | sort | uniq | while read filebase
	do
		sbatch simple_pipeline_PE_hg19_mm9.sbatch $filebase
	done
