source("different_functions.R")  
path1 <- "tsunami_Dohmen_2025_EQ.docx"
#path1 <- "earthquake_table_China_4_E.docx"
# path2 <- "earthquake_table_China_4_N.docx"
#path3 <- "earthquake_table_China_4_SW.docx"
#path4 <- "earthquake_table_China_4_NW.docx"

x1_event <- read_doc(path1)
#x2_event <- read_doc(path2)
#x3_event <- read_doc(path3)
#x4_event <- read_doc(path4)
x_event <- rbind(x1_event)
# x_event <- rbind(x1_event, x2_event, x3_event, x4_event)

chain <- data.frame(data = x_event$Proposed.Encoding)

trigger <- c()
target <- c()
n <- dim(chain)[1]
x <- 0
for (i in 1:n){
  m <- length(strsplit(chain$data[i],split = ";" )[[1]])
  x <- x+m
  for(j in 1:m){
    data_get <-  strsplit(chain$data[i],split = ";" )[[1]][j]
    trigger <- append(trigger,strsplit(data_get,split = "," )[[1]][1])
    target <- append(target,strsplit(data_get,split = "," )[[1]][2])
  }
}

tmp_data <- data.frame(from = trigger,to = target)
library(dplyr)
self_loops <- tmp_data %>% filter(from == to)
tmp_data <- bind_rows(tmp_data, self_loops)  # 确保自环不会被删除
final <- tmp_data %>% distinct(.keep_all = TRUE)


library(igraph)
g <- graph_from_data_frame(final, directed=TRUE)

c <-as_adjacency_matrix(g,sparse=FALSE)

graph <- g
V(graph)$degree <- degree(graph)
V(graph)$degreeIN <- degree(graph, mode = "in")
V(graph)$degreeOUT <- degree(graph, mode = "out") 
V(graph)$closeness <- centr_clo(graph)$res
V(graph)$closenessIN <- centr_clo(graph, mode = "in")$res
V(graph)$closenessOUT <- centr_clo(graph, mode = "out")$res
V(graph)$betweenness <- centr_betw(graph)$res



# 获取节点和边数量
num_nodes <- vcount(graph)
num_edges <- ecount(graph)
degreeIN <- degree(graph, mode = "in")
degreeOUT <- degree(graph, mode = "out") 
#graph_density <- edge_density(graph,loops = FALSE)

# 输出结果
print(paste("节点数量:", num_nodes))
print(paste("边数量:", num_edges))
print(paste(degreeIN))
print(paste(degreeOUT))
print(V(graph))
#print(paste(betweenness))

pdf(paste0("fig_graph_all_diasaters.pdf"),width = 11, height = 10)

set.seed(123)  # 设置随机种子确保可重复性
# par(mfrow=c(2,2), mar=c(0,0,0,0))
layout_saved <- layout_with_fr(graph)
# pal = colorRampPalette(c('blue','green','red'))
# pal = colorRampPalette(c( 'white','cornflowerblue', 'darkseagreen','lightgoldenrod', 'tan1', 
#                           'tomato', 'red', 'orangered3', 'firebrick', 'darkred'))
# 颜色梯度（沿用你的基准红）
pal <- colorRampPalette(c("#FFF886", "#d21d13"))

# --- 可复用绘图函数 -----------------------------------------
plot_centrality_graph <- function(graph, layout, metric_raw,
                                  size_metric  = NULL,   # 用于映射节点大小的指标，默认用 metric_raw
                                  zero_white   = FALSE,  # 指标为 0 的节点是否涂白
                                  label_top_n  = 8,      # 只给最重要的前 N 个节点加标签
                                  size_range   = c(4, 18)) {
  
  # 1) 颜色：按指标归一化映射到色标
  metric <- metric_raw / max(metric_raw)
  graphCol <- pal(500)[as.numeric(cut(metric, breaks = 500))]
  if (zero_white) graphCol[metric_raw == 0] <- "white"
  
  # 2) 大小：让重要节点显著变大，建立层次（sqrt 压缩极端值）
  if (is.null(size_metric)) size_metric <- metric_raw
  s <- sqrt(size_metric)
  if (max(s) == 0) {
    vsize <- rep(mean(size_range), vcount(graph))
  } else {
    vsize <- size_range[1] + (s / max(s)) * (size_range[2] - size_range[1])
  }
  
  # 3) 标签：只显示最重要的前 N 个节点，其余留空，避免文字堆叠
  ord    <- order(metric_raw, decreasing = TRUE)
  keep   <- head(ord, label_top_n)
  vlabel <- rep(NA, vcount(graph))
  vlabel[keep] <- V(graph)$name[keep]
  
  plot(graph, layout = layout,
       # —— 边：弱化视觉权重，让节点成为焦点 ——
       edge.arrow.size = 0.35,
       edge.width      = 0.5,
       edge.curved     = 0.15,
       edge.color      = adjustcolor("grey60", alpha.f = 0.45),  # 半透明，减少叠加糊成一片
       # —— 节点 ——
       vertex.size        = vsize,
       vertex.color       = graphCol,
       vertex.frame.color = adjustcolor("grey30", alpha.f = 0.5),
       # —— 标签：只标重点 ——
       vertex.label        = vlabel,
       vertex.label.family = "Helvetica",
       vertex.label.color  = "black",
       vertex.label.cex    = 1.2,
       vertex.label.font   = 2,        # 加粗
       vertex.label.dist   = 0)
}

# --- 布局：增大斥力 + 多次迭代，从源头上减少节点重叠 ---------
set.seed(123)
layout_saved <- layout_with_fr(graph, niter = 3000)
# 备选：连通度高、簇明显时 KK 布局更舒展
# layout_saved <- layout_with_kk(graph)

# ======== 矢量输出（SVG，可在 Illustrator/Inkscape 编辑） ========
# SVG 不支持多页，故两张图各存一个文件
library(svglite)

# in-degrees：节点大小 = 入度，颜色 = 入度
svglite("fig_graph_indegree.svg", width = 11, height = 10)
plot_centrality_graph(graph, layout_saved,
                      metric_raw = V(graph)$degreeIN,
                      label_top_n = 41)
dev.off()

# out-degrees：节点大小 = 出度，颜色 = 出度，0 出度涂白
svglite("fig_graph_outdegree.svg", width = 11, height = 10)
plot_centrality_graph(graph, layout_saved,
                      metric_raw = V(graph)$degreeOUT,
                      zero_white = TRUE,
                      label_top_n = 41)
dev.off()


pdf("fig_graph_all_disasters.pdf", width = 11, height = 10)

# ======== in-degrees：节点大小 = 入度，颜色 = 入度 ========
plot_centrality_graph(graph, layout_saved,
                      metric_raw = V(graph)$degreeIN,
                      label_top_n = 41)

# ======== out-degrees：节点大小 = 出度，颜色 = 出度，0 出度涂白 ========
plot_centrality_graph(graph, layout_saved,
                      metric_raw = V(graph)$degreeOUT,
                      zero_white = TRUE,
                      label_top_n = 41)

dev.off()

# jpeg("plot.jpg", width = 800, height = 600)
# # par(mfrow=c(2,2), mar=c(0,0,0,0))
# layout_saved <- layout_with_fr(graph)
# #pal = colorRampPalette(c('blue','green','red'))
# pal = colorRampPalette(c('yellow','red'))
# 
# metric <- V(graph)$degreeIN/max(V(graph)$degreeIN)
# graphCol = pal(500)[as.numeric(cut(metric, breaks = 500))]
# plot(graph, layout=layout_saved,
#      edge.arrow.size = .5, edge.curved=.3,
#      vertex.label.family = "Helvetica", vertex.label.color = "black",
#      vertex.color = graphCol,
#      vertex.label.dist = 0, vertex.label.degree = pi/2)
# dev.off()

nodes = num_nodes  

for(i in 1:nodes){
  for(j in 1:nodes){
    if(c[i,j]==1){
      c[i,j] <- 0.1
    }
  }
}
label <- c('0','1','2','3','4','5','6','7','8','9')
# colours <- c('white', 'cornflowerblue', 'darkseagreen',
#              'lightgoldenrod', 'tan1', 'tomato', 'red', 'orangered3', 'firebrick', 'darkred')
# colours <- c('#FFFFFF', '#FFF8E7', '#FFEEB9',
#              '#FFE28A', '#FFD15C', '#FFB72A', '#FF8C00', '#E35502', '#B51A00', '#8B0000')

colours <- c('#FFFFFF', '#FFEEB9', '#FFE28A', 
             '#FFD15C', '#FFB72A', '#FF8C00', '#E35502', '#B51A00', '#9F0000','#8B0000')

# colours <- pal()

pi <- .1

list.ev <- c()
for( i in 1:nodes){
  list.ev <- append(list.ev,colnames(c)[i])
}
for(tau in 1:11){
  M <- array(rowSums(sapply(1:tau, function(i) c %^% i)), dim=c(nodes,nodes))
  rownames(M) <- colnames(M) <- list.ev
  
  A.gg <- melt(c, value.name = "a_ij", varnames = c("Trigger", "Target"))
  A.gg$PrRange <- A.gg$a_ij
  M.gg <- melt(M, value.name = "a_ij", varnames = c("Trigger", "Target"))
  
  A.gg <- mk_colors(A.gg); M.gg <- mk_colors(M.gg)
  
  pdf(paste0('EQ_AND_OTHER_CHINA_tau', tau,'.pdf'))
  gA <- ggplot(A.gg) +
    geom_tile(aes(x = Target, y = factor(Trigger, levels = rev(list.ev)), fill = PrRange), colour = "grey") +
    scale_fill_manual("PrRange", values = colours, limits = label) +
    coord_fixed() +
    theme_minimal() +
    scale_x_discrete(position = "top") +
    theme(axis.text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), legend.position="none") +
    # labs(title = "A", x = "Target j", y = "Trigger i")
    labs( title="A" ,x = "Target j", y = "Trigger i")+
    theme(
      plot.title = element_text(size = 10,face = "bold"),
      plot.title.position = "plot",
      axis.title.x = element_text(size = 8), # 设置x轴标题字体大小
      axis.title.y = element_text(size = 8)  # 设置y轴标题字体大小
    )
  
  gM <- ggplot(M.gg) +
    geom_tile(aes(x = Target, y = factor(Trigger, levels = rev(list.ev)), fill = PrRange), colour = "grey") +
    scale_fill_manual("PrRange", values = colours, limits = label) +
    coord_fixed() +
    theme_minimal() +
    scale_x_discrete(position = "top") +
    theme(axis.text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), legend.position="none") +
    # legend.text = element_text(size = 12) + # 设置图例文本大小为12
    # legend.title = element_text(size = 15) +  # 设置图例标题大小为15
    # guides(fill = guide_legend(override.aes = list(size = 2))) +  # 设置图例键大小为5
    # theme( legend.key.size = unit(55," pt ")) +
    # labs(title = "M", x = "Target j", y = "Trigger i")
    labs(title="B" , x = "Target j", y = "Trigger i") +
    theme(
      plot.title = element_text(size = 10,face = "bold"),
      plot.title.position = "plot",
      axis.title.x = element_text(size = 8), # 设置x轴标题字体大小
      axis.title.y = element_text(size = 8)  # 设置y轴标题字体大小
    )
  
  grid.arrange(gA, gM, ncol=2, nrow=1) 
  dev.off()
}


# M <- array(rowSums(sapply(1:tau, function(i) c %^% 1)), dim=c(nodes,nodes))
# rownames(M) <- colnames(M) <- list.ev
# 
# A.gg <- melt(c, value.name = "a_ij", varnames = c("Trigger", "Target"))
# A.gg$PrRange <- A.gg$a_ij
# M.gg <- melt(M, value.name = "a_ij", varnames = c("Trigger", "Target"))
# 
# A.gg <- mk_colors(A.gg); M.gg <- mk_colors(M.gg)
# 
# pdf(paste0('EQ_AND_OTHER_CHINA_tau', tau,'.pdf'))
# ggplot(A.gg) +
#   geom_tile(aes(x = Target, y = factor(Trigger, levels = rev(list.ev)), fill = PrRange), colour = "grey") +
#   scale_fill_manual("PrRange", values = colours, limits = label) +
#   coord_fixed() +
#   theme_minimal() +
#   scale_x_discrete(position = "top") +
#   theme(axis.text.x = element_text(angle = 90)) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), legend.position="none") +
#   labs(title = "A", x = "Target j", y = "Trigger i")
# gA

