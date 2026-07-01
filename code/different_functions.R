library(docxtractr)
library(ggplot2)
library(reshape2)
library(expm)
library(gridExtra) 
library(Matrix)

pi <- 0.1
read_doc <- function(path){
  doc <-  read_docx(path) # 告诉R数据源的位置并读入内存
  x <- data.frame(docx_extract_tbl(doc, tbl_number = 1))# 提取.docx文档的第1个表格
  return(x)
}
mk_colors <- function(df){
  df$PrRange[df$a_ij == 0] <- 0
  df$PrRange[df$a_ij > 0 & df$a_ij < pi] <- 1
  df$PrRange[df$a_ij > pi-1e-6 & df$a_ij < pi+1e-6] <- 2
  df$PrRange[df$a_ij > pi & df$a_ij <= pi*1.1] <- 3
  df$PrRange[df$a_ij > pi*1.1 & df$a_ij <= pi*1.2] <- 4
  df$PrRange[df$a_ij > pi*1.2 & df$a_ij <= pi*1.3] <- 5
  df$PrRange[df$a_ij > pi*1.3 & df$a_ij <= pi*1.4] <- 6
  df$PrRange[df$a_ij > pi*1.4 & df$a_ij <= pi*1.5] <- 7
  df$PrRange[df$a_ij > pi*1.5 & df$a_ij <= pi*1.6] <- 8
  df$PrRange[df$a_ij > pi*1.6] <- 9
  df$PrRange <- factor(df$PrRange, levels = c(0,seq(9)))
  return(df)
}