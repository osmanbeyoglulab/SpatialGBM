setwd("~/Desktop/Pitt/Project_zhang_gbm_v2/")
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggpointdensity)
library(ggpubr) 
library(ggrepel)
library(matrixStats)
library(patchwork)
library(RColorBrewer)
library(stringr)
library(tidyverse)
library(viridis)
source('aux_correlation.R')

color_palette = c(brewer.pal(n=9, name = "Blues")[1:8] %>% rev(),
                  brewer.pal(n=9, name = "Oranges")[1:8])

rename_rownames <- function(df){
  new_rownames <- sapply(rownames(df), function(x) {
    gsub("_", "-", x)
  })
  rownames(df) <- new_rownames
  return(df)
}

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Bivariate: TF x TF
# ========== ========== ========== ========== ========== ========== ========== ==========

bivariate_tfa_tfa_ravi = read.csv('dataframes/morani_ravi/bivariate_tfa_tfa_per_sample.csv', row.names=1) %>% 
  rename_rownames()

bivariate_tfa_tfa_ren = merge(
  read.csv('dataframes/morani_ren/bivariate_tfa_tfa_per_sample.csv', row.names=1),
  read.csv('dataframes/morani_sonpatki/bivariate_tfa_tfa_per_sample.csv', row.names=1),
  by = "row.names", 
  all = TRUE
) %>% column_to_rownames("Row.names") %>% 
  rename_rownames()

# Find common TF-pathway pairs
common_tfs <- intersect(rownames(bivariate_tfa_tfa_ravi), rownames(bivariate_tfa_tfa_ren))
bivariate_tfa_tfa_ravi_fil = bivariate_tfa_tfa_ravi[common_tfs, ]
bivariate_tfa_tfa_ren_fil = bivariate_tfa_tfa_ren[common_tfs, ]

# Replace NAs with 0
bivariate_tfa_tfa_ravi_fil[is.na(bivariate_tfa_tfa_ravi_fil)] <- 0
bivariate_tfa_tfa_ren_fil[is.na(bivariate_tfa_tfa_ren_fil)] <- 0

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Dotplot / Density plot
# ========== ========== ========== ========== ========== ========== ========== ==========

row_medians_ravi <- rowMedians(
  as.matrix(bivariate_tfa_tfa_ravi_fil[, sapply(bivariate_tfa_tfa_ravi_fil, is.numeric)]), na.rm = TRUE)
row_medians_ren <- rowMedians(
  as.matrix(bivariate_tfa_tfa_ren_fil[, sapply(bivariate_tfa_tfa_ren_fil, is.numeric)]), na.rm = TRUE)

df_moran <- data.frame(
  ravi = row_medians_ravi,
  ren = row_medians_ren,
  row_name = rownames(bivariate_tfa_tfa_ravi_fil)
) 

df_moran_clean <- df_moran %>%
  filter(ravi != 0 & ren != 0)
s = cor.test(df_moran_clean$ravi, df_moran_clean$ren, method = "spearman")
print(s$estimate)

df_moran_dotplot <- df_moran_clean %>%
  mutate(
    label = case_when(
      ravi > 0.55 & ren > 0.55 ~ row_name,
      ravi < -0.6 & ren < -0.6 ~ row_name,
      TRUE ~ ""
    )
  )

p <- ggplot(df_moran_dotplot, aes(x = ravi, y = ren)) +
  geom_abline(intercept = 0, slope = 1, color = "#888888", linetype = "dashed", size = 1) +
  annotate("text", x = 0.5, y = -0.8,
           label = expression(paste("Spearman ", rho, " = ", 0.4305)),
           size = 4, hjust = 0, vjust = 0.5) +
  geom_pointdensity(size = 2) +
  scale_color_viridis(name = "Density") +
  geom_label_repel(
    aes(label = label),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.5,
    label.padding = 0.1,
    label.size = 0.2,
    segment.color = "#666666",
    segment.size = 0.5,
    min.segment.length = 0,
    force = 50,
    force_pull = 2,
    direction = "both",     # 允许标签向各个方向移动
    nudge_x = 0.02,         # 轻微调整初始位置
    nudge_y = 0.02,
    seed = 123              # 设置随机种子，使结果可重复
  ) +
  xlim(-0.8, 1.2) + ylim(-0.8, 1.2) + 
  labs(x = "Discovery cohort (Median Moran's I)", 
       y = "Validation cohort (Median Moran's I)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = "none"
  )

p

my_save_pdf(p, 'plots_manuscript/dot_bimoran_tf_tf.pdf', width = 8, height = 8)

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Heatmap
# ========== ========== ========== ========== ========== ========== ========== ==========

bivariate_tfa_tfa_ravi_subset = filter_dataframe_by_threshold(
  bivariate_tfa_tfa_ravi_fil, thre_val = 0.4, thre_sam = 18)
common_tfs = rownames(bivariate_tfa_tfa_ravi_subset)
common_tfs

tf_groups <- list(
  "SOX2" = c("EPAS1-SOX2", "HIF1A-SOX2", "KLF3-SOX2", "LHX2-SOX2", "NR3C1-SOX2", "PAX6-SOX2", "PRDM1-SOX2", "SMAD3-SOX2"),
  "PAX6" = c("ETV4-PAX6", "HIF1A-PAX6", "PAX6-PRDM1", "PAX6-SMAD3"),
  "HIF1A" = c("HIF1A-SMAD3", "HIF1A-PRDM1"),
  "SMAD3" = c("NR3C1-SMAD3", "PRDM1-SMAD3"),
  "EPAS1" = c("EPAS1-KLF3", "EPAS1-LHX2", "EPAS1-MBD2"),
  "ETV6" = c("ETV6-KLF10", "ETV6-MEF2C", "ETV6-MZF1"),
  "ZNF281" = c("ATF5-ZNF281", "MAFB-ZNF281", "MAZ-ZNF281"),
  "GTF2B" = c("GTF2B-JUND", "GTF2B-KLF3"),
  "Other" = c("ATF5-MAZ", "FOXM1-KLF10", "MEF2C-ZNF711", "MZF1-ZFHX3")
)

desired_order <- unlist(tf_groups)
tf_group_df <- data.frame(
  tf_gene = desired_order,
  cell_type = factor(rep(names(tf_groups), times = sapply(tf_groups, length)),
                     levels = names(tf_groups)))

tfa_long <- bivariate_tfa_tfa_ravi_subset[desired_order, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("tf_gene") %>%
  pivot_longer(-tf_gene, names_to = "sample", values_to = "moran_i") %>%
  left_join(tf_group_df, by = "tf_gene") %>%
  mutate(tf_gene = factor(tf_gene, levels = desired_order))

p_tfa <- ggplot(tfa_long, aes(x = sample, y = tf_gene)) +
  geom_tile(aes(fill = moran_i)) +
  facet_grid(cell_type ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_classic() +
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(
    name = "Bivariate Moran's I",
    colors = color_palette,
    limits = c(-1, 1),
    na.value = "white"
  ) +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background.y = element_rect(fill = "grey95", color = "grey"),
    strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 1),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    text = element_text(size = 10),
    panel.spacing = unit(0.1, "lines"),
    panel.border = element_rect(fill = NA, color = "grey80")
  )

p_tfa

my_save_pdf(p_tfa, 'plots_manuscript/ravi_bimoran_tfa_tfa.pdf', width = 5, height = 6)

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Network
# ========== ========== ========== ========== ========== ========== ========== ==========

library(igraph)

row_medians_ravi <- rowMedians(as.matrix(bivariate_tfa_tfa_ravi_subset[, sapply(bivariate_tfa_tfa_ravi_subset, is.numeric)]),
                               na.rm = TRUE)

edges <- data.frame(
  from = sapply(strsplit(names(row_medians_ravi), "-"), "[", 1),
  to = sapply(strsplit(names(row_medians_ravi), "-"), "[", 2),
  weight = as.numeric(row_medians_ravi)
)

g <- graph_from_data_frame(edges, directed = FALSE)

E(g)$width <- 1 + (E(g)$weight - min(E(g)$weight)) / 
  (max(E(g)$weight) - min(E(g)$weight)) * 4

pdf("plots_manuscript/ravi_network_bimoran_tfa_tfa.pdf", width = 4.5, height = 9)

par(mar = c(2, 2, 2, 2))

# set.seed(123)
# layout_matrix <- layout_with_fr(g, 
#                                 niter = 2000,
#                                 repulserad = 4,
#                                 area = 4000,
#                                 coolexp = 0.95)

plot(g, 
     vertex.size = 25,
     vertex.label.cex = 0.75,
     vertex.label.color = "black",
     vertex.label.dist = 0,
     vertex.label.degree = 0,
     vertex.color = "#A1D99B",
     edge.width = E(g)$width,
     edge.color = "#00441B",
     # layout = layout_matrix,
     main = "",
     frame = TRUE)

legend("bottomleft", 
       legend = round(seq(min(E(g)$weight), max(E(g)$weight), length.out = 4), 2),
       lwd = seq(1, 5, length.out = 4),
       title = "Bivariate\nMoran's I",
       cex = 0.7,
       bty = "n")

dev.off()

brewer.pal(9,'Greens');
"#F7FCF5" "#E5F5E0" "#C7E9C0" "#A1D99B" "#74C476" "#41AB5D" "#238B45" "#006D2C" "#00441B"


components <- decompose(g)

if (5==5){
  # 设置总高度，第二个panel占一半
  pdf("plots_manuscript/ravi_network_bimoran_tfa_tfa.pdf", width = 6.4, height = 8)
  
  # 使用layout矩阵来控制每个panel的高度比例
  # 方法1: 使用layout函数设置高度比例
  layout(matrix(c(1, 2, 3), ncol = 1), heights = c(0.2, 0.55, 0.25))
  
  par(mar = c(1, 1, 1, 1))
  
  # 获取所有components的全局weight范围
  all_weights <- unlist(lapply(components, function(x) E(x)$weight))
  global_min <- min(all_weights)
  global_max <- max(all_weights)
  
  # 统一的legend设置
  legend_breaks <- round(seq(global_min, global_max, length.out = 4), 2)
  legend_lwd <- seq(1, 5, length.out = 4)
  
  # 绘制第一个panel (小)
  g_component <- components[[1]]
  if (global_min != global_max) {
    E(g_component)$width <- 1 + (E(g_component)$weight - global_min) / 
      (global_max - global_min) * 4
  } else {
    E(g_component)$width <- rep(3, length(E(g_component)$weight))
  }
  set.seed(123 + 1)
  layout_matrix <- layout_with_fr(g_component, niter = 2000)
  plot(g_component, 
       vertex.size = 5,
       vertex.label.cex = 1.5,
       vertex.label.color = "black",
       vertex.label.dist = 0,
       vertex.label.degree = 0,
       vertex.color = "#E5F5E0",
       edge.width = E(g_component)$width,
       edge.color = "#A1D99B",
       layout = layout_matrix,
       frame = TRUE)
  
  # 绘制第二个panel (占一半高度)
  g_component <- components[[2]]
  if (global_min != global_max) {
    E(g_component)$width <- 1 + (E(g_component)$weight - global_min) / 
      (global_max - global_min) * 4
  } else {
    E(g_component)$width <- rep(3, length(E(g_component)$weight))
  }
  set.seed(123 + 2)
  layout_matrix <- layout_with_fr(g_component, niter = 2000)
  plot(g_component, 
       vertex.size = 1.25,
       vertex.label.cex = 1.5,
       vertex.label.color = "black",
       vertex.label.dist = 0,
       vertex.label.degree = 0,
       vertex.color = "#E5F5E0",
       edge.width = E(g_component)$width,
       edge.color = "#A1D99B",
       layout = layout_matrix,
       frame = TRUE)
  
  # 绘制第三个panel (小)
  g_component <- components[[3]]
  if (global_min != global_max) {
    E(g_component)$width <- 1 + (E(g_component)$weight - global_min) / 
      (global_max - global_min) * 4
  } else {
    E(g_component)$width <- rep(3, length(E(g_component)$weight))
  }
  set.seed(123 + 3)
  layout_matrix <- layout_with_fr(g_component, niter = 2000)
  plot(g_component, 
       vertex.size = 5,
       vertex.label.cex = 1.5,
       vertex.label.color = "black",
       vertex.label.dist = 0,
       vertex.label.degree = 0,
       vertex.color = "#E5F5E0",
       edge.width = E(g_component)$width,
       edge.color = "#A1D99B",
       layout = layout_matrix,
       frame = TRUE)
  
  # 在第三个panel添加legend
  legend("bottomleft", 
         legend = legend_breaks,
         lwd = legend_lwd,
         title = "Bivariate\nMoran's I",
         cex = 1,
         bty = "n")
  
  dev.off()
}

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Validation
# ========== ========== ========== ========== ========== ========== ========== ==========

bivariate_tfa_tfa_ren_subset = filter_dataframe_by_threshold(
  bivariate_tfa_tfa_ren_fil, thre_val = 0.3, thre_sam = 13)
common_tfs_ren = rownames(bivariate_tfa_tfa_ren_subset) %>% intersect(rownames(bivariate_tfa_tfa_ravi_subset))
common_tfs_ren

tf_groups_ren <- list(
  "SOX2" = c("LHX2-SOX2", "PRDM1-SOX2", "SMAD3-SOX2"),
  "PAX6" = c("HIF1A-PAX6", "PAX6-PRDM1", "PAX6-SMAD3"), 
  "SMAD3" = c("PRDM1-SMAD3", "HIF1A-SMAD3"),  
  "Other" = c("MEF2C-ZNF711", "MAZ-ZNF281") 
)

desired_order_ren <- unlist(tf_groups_ren)
tf_group_df_ren <- data.frame(
  tf_gene = desired_order_ren,
  cell_type = factor(rep(names(tf_groups_ren), times = sapply(tf_groups_ren, length)),
                     levels = names(tf_groups_ren)))

tfa_long_ren <- bivariate_tfa_tfa_ren_subset[desired_order_ren, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("tf_gene") %>%
  pivot_longer(-tf_gene, names_to = "sample", values_to = "moran_i") %>%
  left_join(tf_group_df_ren, by = "tf_gene") %>%
  mutate(tf_gene = factor(tf_gene, levels = desired_order_ren))

p_tfa_ren <- ggplot(tfa_long_ren, aes(x = sample, y = tf_gene)) +
  geom_tile(aes(fill = moran_i)) +
  facet_grid(cell_type ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_classic() +
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(
    name = "Bivariate Moran's I",
    colors = color_palette,
    limits = c(-1, 1),
    na.value = "white"
  ) +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background.y = element_rect(fill = "grey95", color = "grey"),
    strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 1),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    text = element_text(size = 10),
    panel.spacing = unit(0.1, "lines"),
    panel.border = element_rect(fill = NA, color = "grey80")
  )

p_tfa_ren

my_save_pdf(p_tfa_ren, 'plots_manuscript/ren_bimoran_tfa_tfa.pdf', width = 4.8, height = 3.6)

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Network
# ========== ========== ========== ========== ========== ========== ========== ==========

bivariate_tfa_tfa_ren_subset = bivariate_tfa_tfa_ren_subset[common_tfs_ren, ]
row_medians_ren <- rowMedians(as.matrix(bivariate_tfa_tfa_ren_subset[, sapply(bivariate_tfa_tfa_ren_subset, is.numeric)]),
                               na.rm = TRUE)

edges <- data.frame(
  from = sapply(strsplit(names(row_medians_ren), "-"), "[", 1),
  to = sapply(strsplit(names(row_medians_ren), "-"), "[", 2),
  weight = as.numeric(row_medians_ren)
)

g <- graph_from_data_frame(edges, directed = FALSE)

E(g)$width <- 1 + (E(g)$weight - min(E(g)$weight)) / 
  (max(E(g)$weight) - min(E(g)$weight)) * 4

pdf("plots_manuscript/ren_network_bimoran_tfa_tfa.pdf", width = 4.8, height = 3.6)

par(mar = c(2, 2, 2, 2))

# set.seed(123)
# layout_matrix <- layout_with_fr(g, 
#                                 niter = 2000,
#                                 repulserad = 4,
#                                 area = 4000,
#                                 coolexp = 0.95)

plot(g, 
     vertex.size = 15,
     vertex.label.cex = 0.75,
     vertex.label.color = "black",
     vertex.label.dist = 0,
     vertex.label.degree = 0,
     vertex.color = "#A1D99B",
     edge.width = E(g)$width,
     edge.color = "#00441B",
     # layout = layout_matrix,
     main = "",
     frame = TRUE)

legend("bottomleft", 
       legend = round(seq(min(E(g)$weight), max(E(g)$weight), length.out = 4), 2),
       lwd = seq(1, 5, length.out = 4),
       title = "Bivariate\nMoran's I",
       cex = 0.7,
       bty = "n")

dev.off()
