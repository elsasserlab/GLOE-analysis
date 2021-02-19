library(elsasserlib)
library(ggplot2)
library(reshape2)
library(ggpubr)

setwd("~/Downloads/GLOE/210121/")

samples <- c("WT-UT-1-1","WT-UT-1-2","WT-UT-1-3","WT-20-2-1","WT-20-2-2","WT-20-2-3","WT-24hr-3-1","WT-24hr-3-2","WT-24hr-3-3","KO-UT-4-1","KO-UT-4-2","KO-UT-4-3","KO-20-5-1","KO-20-5-2","KO-20-5-3","KO-24hr-6-1","KO-24hr-6-2","KO-24hr-6-3")

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

plot_bw_bed_summary_heatmap(fn$files,"HCT116_ChromHMM.hg19.chr22.bed",labels = fn$labels)

plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_ATACseq.hg19.bed")
plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_CTCF.hg19.bed")
plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_KAP1.hg19.bed")
plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_MYC.hg19.bed")
plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_TOPI.hg19.bed")
plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_YY1.hg19.bed")
plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_RNAP2.hg19.bed")

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

plot_bw_bins_violin(bwfiles = fn$files,labels = fn$labels,bin_size = 10000)

plot_bw_bins_scatter(x=fn$files[1],y=fn$files[2],bin_size = 10000)

plot_bw_profile(bwfiles = fn$files,labels = fn$labels,mode = "center",upstream = 1000, downstream = 1000,bedfile = "HCT116_ATACseq.hg19.bed")
