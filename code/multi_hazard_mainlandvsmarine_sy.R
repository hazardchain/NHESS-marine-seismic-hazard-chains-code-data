## ============================================================
## 海洋地震链 vs 大陆地震链：建图 + 单网绘图 + 统一 15 类对比
## ============================================================
source("different_functions.R")
library(dplyr)
library(igraph)
library(svglite)

## ============================================================
## 0. 通用函数
## ============================================================

## 0a. 从一组 docx 构建灾害链有向图（封装原 pipeline）
build_graph <- function(paths){
  x_event <- do.call(rbind, lapply(paths, read_doc))
  chain   <- data.frame(data = x_event$Proposed.Encoding)
  
  trigger <- c(); target <- c()
  for (i in seq_len(nrow(chain))){
    parts <- strsplit(chain$data[i], split = ";")[[1]]
    for (p in parts){
      pair <- strsplit(p, split = ",")[[1]]
      trigger <- append(trigger, trimws(pair[1]))   # trimws 防止 " TS" 与 "TS" 被当成两类
      target  <- append(target,  trimws(pair[2]))
    }
  }
  tmp        <- data.frame(from = trigger, to = target)
  self_loops <- tmp %>% filter(from == to)
  tmp        <- bind_rows(tmp, self_loops)           # 保留自环
  final      <- tmp %>% distinct(.keep_all = TRUE)
  
  g <- graph_from_data_frame(final, directed = TRUE)
  V(g)$degreeIN  <- degree(g, mode = "in")
  V(g)$degreeOUT <- degree(g, mode = "out")
  g
}

## 0b. 单网度中心性绘图函数（黄→红配色，标签居中）
pal <- colorRampPalette(c("#FFF886", "#d21d13"))
plot_centrality_graph <- function(graph, layout, metric_raw,
                                  size_metric = NULL, zero_white = FALSE,
                                  label_top_n = Inf, size_range = c(4, 18),
                                  label_cex = 1.2){
  metric   <- metric_raw / max(metric_raw)
  graphCol <- pal(500)[as.numeric(cut(metric, breaks = 500))]
  if (zero_white) graphCol[metric_raw == 0] <- "white"
  
  if (is.null(size_metric)) size_metric <- metric_raw
  s <- sqrt(size_metric)
  vsize <- if (max(s) == 0) rep(mean(size_range), vcount(graph))
  else size_range[1] + (s / max(s)) * (size_range[2] - size_range[1])
  
  keep   <- head(order(metric_raw, decreasing = TRUE), label_top_n)
  vlabel <- rep(NA, vcount(graph)); vlabel[keep] <- V(graph)$name[keep]
  
  plot(graph, layout = layout,
       edge.arrow.size = 0.35, edge.width = 0.5, edge.curved = 0.15,
       edge.color = adjustcolor("grey60", alpha.f = 0.45),
       vertex.size = vsize, vertex.color = graphCol,
       vertex.frame.color = adjustcolor("grey30", alpha.f = 0.5),
       vertex.label = vlabel, vertex.label.family = "Helvetica",
       vertex.label.color = "black", vertex.label.cex = label_cex,
       vertex.label.font = 2, vertex.label.dist = 0)
}

## 0c. 单网络出图（in / out 两张，SVG + PDF）
draw_network <- function(graph, prefix, w = 11, h = 10){
  set.seed(123)
  lay <- layout_with_fr(graph, niter = 3000)
  
  svglite(paste0(prefix, "_indegree.svg"), width = w, height = h)
  plot_centrality_graph(graph, lay, V(graph)$degreeIN)
  dev.off()
  svglite(paste0(prefix, "_outdegree.svg"), width = w, height = h)
  plot_centrality_graph(graph, lay, V(graph)$degreeOUT, zero_white = TRUE)
  dev.off()
  
  pdf(paste0(prefix, ".pdf"), width = w, height = h)
  plot_centrality_graph(graph, lay, V(graph)$degreeIN)
  plot_centrality_graph(graph, lay, V(graph)$degreeOUT, zero_white = TRUE)
  dev.off()
}

## ============================================================
## 1. 构建两张网络
## ============================================================
marine_paths   <- "tsunami_Dohmen_2025_EQ.docx"        # ← 确认海洋链文件名
mainland_paths <- c("earthquake_table_China_4_E.docx",
                    "earthquake_table_China_4_N.docx",
                    "earthquake_table_China_4_SW.docx",
                    "earthquake_table_China_4_NW.docx")

g_marine   <- build_graph(marine_paths)
g_mainland <- build_graph(mainland_paths)

for (nm in c("g_marine", "g_mainland")){
  g <- get(nm)
  cat(sprintf("\n[%s] 节点=%d  边=%d\n", nm, vcount(g), ecount(g)))
  print(V(g)$name)
}

## ============================================================
## 2. 各自单网绘图（原 pipeline 的图，每网 in/out 两张）
## ============================================================
draw_network(g_marine,   "fig_marine")
draw_network(g_mainland, "fig_mainland")

## ============================================================
## 3. 海洋链 39 → 15 类合并（与大陆链口径对齐）
## ============================================================
coarse_map <- c(
  SEQ="EQ or AS", AS="EQ or AS", EQ="EQ or AS", VE="EQ or AS", AI="EQ or AS",
  GU="Mass Slide", SLS="Mass Slide", CLS="Mass Slide", LI="Mass Slide",
  MB="Mass Slide", SUB="Mass Slide", LS="Mass Slide", RO="Mass Slide",
  LB="Mass Slide", SA="Mass Slide",
  PHF="Critical Infra. Fail.", SDF="Critical Infra. Fail.",
  CF="Critical Infra. Fail.", NFA="Critical Infra. Fail.", SPF="Critical Infra. Fail.",
  BI="Business Interr.", MTD="Business Interr.",
  WP="Water Dil.", WD="Water Dil.",
  NF="Critical NET Fail.", SCF="Critical NET Fail.",
  FA="Flora an Fauna", MED="Flora an Fauna",
  FI="Fire", EX="Fire",
  IR="Healthcare Deg.", HD="Healthcare Deg.",
  # Flood = 广义水文/气象淹没：海啸、沿海洪涝、河洪、风暴潮、台风及其他风暴
  TS="Flood", CFL="Flood", FF="Flood", SS="Flood", TC="Flood",
  WS="Flood", RS="Flood", OS="Flood", EW="Flood",
  EC="Econ. Crisis", DI="Disease", PP="Phy.-Psych. Trauma",
  CP="Culture Perils", SU="Social Unrest")

miss <- setdiff(V(g_marine)$name, names(coarse_map))
if (length(miss)) stop("以下海洋链节点未在 coarse_map 中映射：", paste(miss, collapse=", "))

grp <- factor(coarse_map[V(g_marine)$name])
g_marine_c <- contract(g_marine, as.integer(grp), vertex.attr.comb = "first")
V(g_marine_c)$name <- levels(grp)
## 去重边 + 去自环，使两网口径一致
g_marine_c   <- simplify(g_marine_c,   remove.multiple = TRUE, remove.loops = TRUE)
g_mainland_c <- simplify(g_mainland,   remove.multiple = TRUE, remove.loops = TRUE)

## ============================================================
## 4. 可比拓扑指标
## ============================================================
topo <- function(g, name){
  ind <- degree(g, mode="in"); outd <- degree(g, mode="out")
  data.frame(network=name, nodes=vcount(g), edges=ecount(g),
             density  = round(edge_density(g, loops=FALSE), 3),
             mean_deg = round(ecount(g)/vcount(g), 2),
             max_out=max(outd), max_in=max(ind),
             Cout = round(centr_degree(g, mode="out", loops=FALSE)$centralization, 3),
             Cin  = round(centr_degree(g, mode="in",  loops=FALSE)$centralization, 3),
             sources=sum(ind==0), sinks=sum(outd==0),
             mean_path=round(mean_distance(g, directed=TRUE), 2))
}
cat("\n===== 可比拓扑指标 =====\n")
print(rbind(topo(g_marine_c, "Marine(coarse)"),
            topo(g_mainland_c, "Mainland")), row.names = FALSE)

## ============================================================
## 5. 四层能量瀑布
## ============================================================
layer15 <- c(
  "EQ or AS"="Trigger",
  "Mass Slide"="Natural", "Water Dil."="Natural", "Flora an Fauna"="Natural",
  "Flood"="Natural", "Disease"="Natural",
  "Critical Infra. Fail."="Engineering", "Critical NET Fail."="Engineering", "Fire"="Engineering",
  "Business Interr."="Consequence", "Healthcare Deg."="Consequence",
  "Econ. Crisis"="Consequence", "Phy.-Psych. Trauma"="Consequence",
  "Culture Perils"="Consequence", "Social Unrest"="Consequence")
L <- c("Trigger", "Natural", "Engineering", "Consequence")

waterfall <- function(g, name){
  lay <- factor(layer15[V(g)$name], levels = L)
  o <- tapply(degree(g, mode="out"), lay, sum); o[is.na(o)] <- 0
  i <- tapply(degree(g, mode="in"),  lay, sum); i[is.na(i)] <- 0
  data.frame(network=name, layer=L, out=as.integer(o[L]), inn=as.integer(i[L]),
             net=as.integer(o[L]-i[L]),
             net_pct=round((o[L]-i[L])/ecount(g)*100, 1))
}
cat("\n===== 海洋链(粗粒度) 层级瀑布 =====\n")
print(waterfall(g_marine_c, "Marine(coarse)"), row.names = FALSE)
cat("\n===== 大陆链 层级瀑布 =====\n")
print(waterfall(g_mainland_c, "Mainland"), row.names = FALSE)

## ============================================================
## 6. 统一 15 类、统一布局的 2×2 对比图
## ============================================================
abbr <- c("EQ or AS"="EQ", "Mass Slide"="MS", "Critical Infra. Fail."="CIF",
          "Business Interr."="BI", "Water Dil."="WD", "Critical NET Fail."="CNF",
          "Flora an Fauna"="FA", "Fire"="FI", "Healthcare Deg."="HD", "Flood"="FL",
          "Econ. Crisis"="EC", "Disease"="DI", "Phy.-Psych. Trauma"="PP",
          "Culture Perils"="CP", "Social Unrest"="SU")

common <- sort(unique(c(V(g_marine_c)$name, V(g_mainland_c)$name)))
uni <- graph_from_data_frame(
  rbind(igraph::as_data_frame(g_marine_c), igraph::as_data_frame(g_mainland_c)),
  directed = TRUE, vertices = data.frame(name = common))
set.seed(123)
lay_u <- layout_with_fr(uni, niter = 3000); rownames(lay_u) <- V(uni)$name

plot_one <- function(g, mode, title){
  m_raw  <- degree(g, mode = mode)
  metric <- if (max(m_raw)==0) rep(0, length(m_raw)) else m_raw/max(m_raw)
  col    <- pal(500)[as.numeric(cut(metric, breaks = 500))]
  if (mode=="out") col[m_raw==0] <- "white"
  s <- sqrt(m_raw); vs <- 12 + (if (max(s)==0) 0 else s/max(s))*26
  plot(g, layout = lay_u[V(g)$name, ],
       vertex.size=vs, vertex.color=col,
       vertex.frame.color=adjustcolor("grey30", alpha.f=0.5),
       vertex.label=abbr[V(g)$name], vertex.label.family="Helvetica",
       vertex.label.color="black", vertex.label.cex=1.5,
       vertex.label.font=2, vertex.label.dist=0,
       edge.color=adjustcolor("grey60", alpha.f=0.5), edge.width=0.7,
       edge.arrow.size=0.35, edge.curved=0.15, main=title)
}
draw_panel <- function(){
  par(mfrow=c(2,2), mar=c(1,1,2.5,1))
  plot_one(g_marine_c,   "in",  "Marine (coarse) — in-degree")
  plot_one(g_marine_c,   "out", "Marine (coarse) — out-degree")
  plot_one(g_mainland_c, "in",  "Mainland — in-degree")
  plot_one(g_mainland_c, "out", "Mainland — out-degree")
}
svglite("fig_compare_coarse.svg", width=14, height=13); draw_panel(); dev.off()
pdf("fig_compare_coarse.pdf",     width=14, height=13); draw_panel(); dev.off()

cat("\n完成。输出文件：\n",
    " 单网图: fig_marine_*.svg/.pdf, fig_mainland_*.svg/.pdf\n",
    " 对比图: fig_compare_coarse.svg/.pdf\n")


## ============================================================
## 7. 汇总为论文 Table 4（Marine vs Mainland 两列）
## ============================================================
build_table4 <- function(g_marine_c, g_mainland_c){
  tm <- topo(g_marine_c,   "Marine")
  tl <- topo(g_mainland_c, "Mainland")
  wm <- waterfall(g_marine_c,   "Marine")
  wl <- waterfall(g_mainland_c, "Mainland")
  
  ## 取各层 net_pct（按 Trigger/Natural/Engineering/Consequence 顺序）
  pm <- setNames(wm$net_pct, wm$layer)
  pl <- setNames(wl$net_pct, wl$layer)
  sgn <- function(x) sprintf("%+.1f", x)     # 带正负号，如 +19.0 / -31.0
  
  metric <- c("Nodes / Edges",
              "Network density",
              "Out-degree centralization",
              "In-degree centralization",
              "Max out-degree",
              "Pure sources / sinks",
              "Mean path length",
              "Net flow — Trigger layer (%)",
              "Net flow — Natural layer (%)",
              "Net flow — Engineering layer (%)",
              "Net flow — Consequence layer (%)")
  
  higher <- c("—",
              "More hazards are interconnected",
              "Triggering concentrated in one source",
              "Consequences converge more sharply",
              "A more dominant single trigger",
              "More terminal endpoints",
              "Almost equal",
              "Stronger source layer",
              "(near balance = conduit)",
              "Forwarding",
              "Sharper terminal absorption")
  
  col <- function(t, p){
    c(paste(t$nodes, "/", t$edges),
      format(t$density,  nsmall = 3),
      format(t$Cout,     nsmall = 3),
      format(t$Cin,      nsmall = 3),
      t$max_out,
      paste(t$sources, "/", t$sinks),
      format(t$mean_path, nsmall = 2),
      sgn(p["Trigger"]), sgn(p["Natural"]),
      sgn(p["Engineering"]), sgn(p["Consequence"]))
  }
  
  data.frame(Metric   = metric,
             Marine   = col(tm, pm),
             Mainland = col(tl, pl),
             `Higher value indicates` = higher,
             check.names = FALSE, row.names = NULL)
}

tab4 <- build_table4(g_marine_c, g_mainland_c)
cat("\n===== Table 4: Marine vs Mainland =====\n")
print(tab4, row.names = FALSE, right = FALSE)
write.csv(tab4, "table4_marine_vs_mainland.csv", row.names = FALSE)
