library(wigglescout)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(GenomicRanges)

setwd("~/Downloads/GLOE/210121/")

replicates <- c("WT-UT-1-1","WT-UT-1-2","WT-UT-1-3","WT-20-2-1","WT-20-2-2","WT-20-2-3","WT-24hr-3-1","WT-24hr-3-2","WT-24hr-3-3","KO-UT-4-1","KO-UT-4-2","KO-UT-4-3","KO-20-5-1","KO-20-5-2","KO-20-5-3","KO-24hr-6-1","KO-24hr-6-2","KO-24hr-6-3")
samples <- c("WT-UT-1","WT-20-2","WT-24hr-3","KO-UT-4","KO-20-5","KO-24hr-6")

build_filenames <- function(path, cond) {
  files <- paste0(path,"/GL-",cond,".hg19.bs50.bw")
  labels <- paste0("GLOE ", cond)
  return(data.frame(files=files,labels=labels))
}

build_filenames_plus <- function(path, cond) {
  files <- paste0(path,"/GL-",cond,".hg19.R1plus.bw")
  labels <- paste0("GLOE (+) ", cond)
  return(data.frame(files=files,labels=labels))
}

build_filenames_minus <- function(path, cond) {
  files <- paste0(path,"//GL-",cond,".hg19.R1minus.bw")
  labels <- paste0("GLOE (-) ", cond)
  return(data.frame(files=files,labels=labels))
}

fn <- build_filenames("bw/",samples)

plot_bw_loci_summary_heatmap(fn$files,"HCT116_WT_10_segments.bed",labels = fn$labels, aggregate_by = "median")


explot_bw_loci_scatter(x=fn$files[1],y=fn$files[6],loci = "1kb.hg19.chr10.bed", highlight = "HCT116_CTCF.hg19.bed") +
  scale_x_continuous(trans = 'log2', limits = c(1,100)) +
  scale_y_continuous(trans = 'log2', limits = c(1,100))
 
plot_bw_loci_scatter(x=fn$files[1],y=fn$files[3],loci = "HCT116_CTCF.hg19.bed",) +
  scale_x_continuous(trans = 'log2', limits = c(1,100)) +
  scale_y_continuous(trans = 'log2', limits = c(1,100))

df <- bw_loci(bwfiles = fn$files, labels = fn$labels,loci = "HCT116_CTCF.hg19.bed",remove_top = 0.001)

df_1 <- select(data.frame(df),contains("GLOE"))
mdf <- melt(df_1)

ggviolin(mdf,"variable","value") + scale_y_continuous(trans = 'log2')
ggboxplot(mdf,"variable","value") + scale_y_continuous(trans = 'log2')

ggpaired(df_1,cond1="GLOE.WT.UT.1",cond2="GLOE.WT.24hr.3", line.color = "gray", line.size = 0.4, palette = "jco") + 
  stat_compare_means(paired = TRUE)

plot_bw_profile(bwfiles = fn$files, labels = fn$labels,bedfile  = "HCT116_CTCF.hg19.bed")

plot_bw_profile(bwfiles = fn$files, labels = fn$labels,bedfile  = "ctcf_CAD_up.bed")

df_1 <- data.frame(df)

CAD <- df_1[(df$GLOE.WT.24hr.3 > max(df$GLOE.WT.UT.1,df$GLOE.WT.20.2, df$GLOE.KO.UT.4, df$GLOE.KO.20.5, df$GLOE.KO.24hr.6)),]

write.table(CAD,file = "ctcf_CAD_up.bed",quote = F, sep = "\t", row.names = F)

fn <- build_filenames_plus("bw/",samples)

plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_ATACseq.hg19.bed")


plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_CTCF.hg19.bed",bin_size = 10)

plot_bw_bed_summary_heatmap(fn$files,"HCT116_ChromHMM.hg19.chr22.bed",labels = fn$labels)


fn <- build_filenames_plus("bw/",samples)

plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_CTCF.hg19.bed",bin_size = 10)

fn <- build_filenames_minus("bw/",samples)

plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_CTCF.hg19.bed",bin_size = 10)
plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_ATACseq.hg19.bed",bin_size = 10)
plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_KAP1.hg19.bed",bin_size = 10)
plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_MYC.hg19.bed",bin_size = 10)
plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_TOPI.hg19.bed",bin_size = 10)
plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_YY1.hg19.bed",bin_size = 10)
plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_RNAP2.hg19.bed",bin_size = 10)




plot_bw_bed_summary_heatmap(fn$files,"HeLa.segHMM.hg19.chr22.bed",labels = fn$labels)

plot_bw_bed_summary_heatmap(fn$files,"Gm12878.segHMM.hg19.chr22.bed",labels = fn$labels)

plot_bw_bed_summary_heatmap(fn$files,"repmasker_sel.hg19.bed",labels = fn$labels)

plot_bw_bed_summary_heatmap(fn$files,"repmasker_cluster_lt1kb.hg19.bed",labels = fn$labels)

plot_bw_bed_summary_heatmap(fn$files,"repmasker.downsamples.hg19.bed",labels = fn$labels)

