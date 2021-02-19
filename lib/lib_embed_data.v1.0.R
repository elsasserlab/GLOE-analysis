if (!require("xfun")) install.packages("xfun")

embed_last_plot=function(name) {
  fn_pdf <- paste0(name,".pdf")
  
  ggsave(filename = fn_pdf, device = "pdf")
  
  xfun::embed_file(fn_pdf)
}

embed_last_data=function(name) {
  fn_tbl <- paste0(name,".txt")
  
  write.table(x = last_plot()$data, file = fn_tbl,quote = F, sep = "\t",row.names = F,col.names = T)
  
  xfun::embed_file(fn_tbl)
}

embed_plot=function(p, name) {
  fn_pdf <- paste0(name,".pdf")
  
  ggsave(filename = fn_pdf, plot = p, device = "pdf")
  
  xfun::embed_file(fn_pdf)
}

embed_data=function(p, name) {
  fn_tbl <- paste0(name,".txt")
  
  write.table(x = p$data, file = fn_tbl,quote = F, sep = "\t",row.names = F,col.names = T)
  
  xfun::embed_file(fn_tbl)
}