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
#           Bivariate: TF x marker
# ========== ========== ========== ========== ========== ========== ========== ==========

bivariate_tfa_marker_ravi = read.csv('dataframes/morani_ravi/bivariate_tfa_marker_per_sample.csv', row.names=1) %>% 
  rename_rownames()

bivariate_tfa_marker_ren = merge(
  read.csv('dataframes/morani_ren/bivariate_tfa_marker_per_sample.csv', row.names=1),
  read.csv('dataframes/morani_sonpatki/bivariate_tfa_marker_per_sample.csv', row.names=1),
  by = "row.names", 
  all = TRUE
) %>% column_to_rownames("Row.names") %>% 
  rename_rownames()

# Find common TF-pathway pairs
common_tfs <- intersect(rownames(bivariate_tfa_marker_ravi), rownames(bivariate_tfa_marker_ren))
bivariate_tfa_marker_ravi_fil = bivariate_tfa_marker_ravi[common_tfs, ]
bivariate_tfa_marker_ren_fil = bivariate_tfa_marker_ren[common_tfs, ]

# Replace NAs with 0
bivariate_tfa_marker_ravi_fil[is.na(bivariate_tfa_marker_ravi_fil)] <- 0
bivariate_tfa_marker_ren_fil[is.na(bivariate_tfa_marker_ren_fil)] <- 0

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Dotplot / Density plot
# ========== ========== ========== ========== ========== ========== ========== ==========

row_medians_ravi <- rowMedians(
  as.matrix(bivariate_tfa_marker_ravi_fil[, sapply(bivariate_tfa_marker_ravi_fil, is.numeric)]), na.rm = TRUE)
row_medians_ren <- rowMedians(
  as.matrix(bivariate_tfa_marker_ren_fil[, sapply(bivariate_tfa_marker_ren_fil, is.numeric)]), na.rm = TRUE)

df_moran <- data.frame(
  ravi = row_medians_ravi,
  ren = row_medians_ren,
  row_name = rownames(bivariate_tfa_marker_ravi_fil)
) 

df_moran_clean <- df_moran %>%
  filter(ravi != 0 & ren != 0)
s = cor.test(df_moran_clean$ravi, df_moran_clean$ren, method = "spearman")
print(s$estimate)

df_moran_dotplot <- df_moran_clean %>%
  mutate(
    label = case_when(
      ravi > 0.4 & ren > 0.4 ~ row_name,
      ravi < -0.4 & ren < -0.4 ~ row_name,
      TRUE ~ ""
    )
  )

p <- ggplot(df_moran_dotplot, aes(x = ravi, y = ren)) +
  geom_abline(intercept = 0, slope = 1, color = "#888888", linetype = "dashed", size = 1) +
  annotate("text", x = 0.4, y = -0.6,
           label = expression(paste("Spearman ", rho, " = ", 0.6062)),
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
  xlim(-0.6, 0.8) + ylim(-0.6, 0.8) + 
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

my_save_pdf(p, 'plots_manuscript/dot_bimoran_tf_mg.pdf', width = 8, height = 8)

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Heatmap
# ========== ========== ========== ========== ========== ========== ========== ==========

bivariate_tfa_marker_ravi_subset = filter_dataframe_by_threshold(
  bivariate_tfa_marker_ravi_fil, thre_val = 0.25, thre_sam = 15)
common_tfs = rownames(bivariate_tfa_marker_ravi_subset)
common_tfs

tfa_long <- bivariate_tfa_marker_ravi_subset %>%
  as.data.frame() %>%
  rownames_to_column("tf_gene") %>%
  pivot_longer(cols = -tf_gene,
               names_to = "sample",
               values_to = "moran_i")
tfa_long$gene <- sapply(strsplit(tfa_long$tf_gene, "-"), function(x) paste(x[-1], collapse = "-"))

gene_to_celltype <- list()
for(cell_type in names(marker_genes)) {
  for(gene in marker_genes[[cell_type]]) {
    gene_to_celltype[[gene]] <- cell_type
  }
}

tfa_long$cell_type <- sapply(tfa_long$gene, function(g) {
  if(g %in% names(gene_to_celltype)) {
    gene_to_celltype[[g]]
  } else {
    "Other"
  }
})

celltype_order <- names(marker_genes)
tfa_long$cell_type <- factor(tfa_long$cell_type,
                                   levels = c(celltype_order, "Other"))

p <- ggplot(tfa_long, aes(x = sample, y = tf_gene)) +
  geom_tile(aes(fill = moran_i)) +  # 改用geom_tile，用moran_i填充
  facet_grid(cell_type ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_classic() +
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(  # 改为scale_fill_gradientn
    name = "TF-Marker\nBivariate\nMoran's I",
    # colors = brewer.pal(n = 11, name = "RdBu") %>% rev(),
    colors = color_palette,
    values = scales::rescale(c(-0.6, -0.3, 0, 0.3, 0.6)),
    limits = c(-0.8, 0.8),
    na.value = "grey"
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
    
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = rel(0.8)),
    legend.text = element_text(size = rel(0.8)),
    text = element_text(size = 12),
    
    # 添加heatmap特定的美化选项
    panel.spacing = unit(0.1, "lines"),  # 减少facet之间的间距
    panel.border = element_rect(fill = NA, color = "grey80"),  # 添加边框区分facet
  )
p
my_save_pdf(p, 'plots_manuscript/ravi_bimoran_tfa_mge.pdf', width = 8, height = 9)

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Validation
# ========== ========== ========== ========== ========== ========== ========== ==========

# bivariate_tfa_marker_ren_subset = filter_dataframe_by_threshold(
#   bivariate_tfa_marker_ren_fil, thre_val = 0.25, thre_sam = 13)
# common_tfs_ren = rownames(bivariate_tfa_marker_ren_subset) %>% intersect(rownames(bivariate_tfa_marker_ravi_subset))
# common_tfs_ren
# 
# bivariate_tfa_marker_ren_subset = filter_dataframe_by_threshold(
#   bivariate_tfa_marker_ren_fil, thre_val = 0.3, thre_sam = 13)
# common_tfs_ren = rownames(bivariate_tfa_marker_ren_subset) 
# common_tfs_ren
# 
# 
# tfa_long <- bivariate_tfa_marker_ravi_subset %>%
#   as.data.frame() %>%
#   rownames_to_column("tf_gene") %>%
#   pivot_longer(cols = -tf_gene,
#                names_to = "sample",
#                values_to = "moran_i")
# tfa_long$gene <- sapply(strsplit(tfa_long$tf_gene, "-"), function(x) paste(x[-1], collapse = "-"))
# 
# gene_to_celltype <- list()
# for(cell_type in names(marker_genes)) {
#   for(gene in marker_genes[[cell_type]]) {
#     gene_to_celltype[[gene]] <- cell_type
#   }
# }
# 
# tfa_long$cell_type <- sapply(tfa_long$gene, function(g) {
#   if(g %in% names(gene_to_celltype)) {
#     gene_to_celltype[[g]]
#   } else {
#     "Other"
#   }
# })
# 
# celltype_order <- names(marker_genes)
# tfa_long$cell_type <- factor(tfa_long$cell_type,
#                              levels = c(celltype_order, "Other"))
# 
# p <- ggplot(tfa_long, aes(x = sample, y = tf_gene)) +
#   geom_tile(aes(fill = moran_i)) +  # 改用geom_tile，用moran_i填充
#   facet_grid(cell_type ~ ., scales = "free_y", space = "free_y", switch = "y") +
#   theme_classic() +
#   scale_x_discrete(position = "top") +
#   scale_fill_gradientn(  # 改为scale_fill_gradientn
#     name = "TF-Marker\nBivariate\nMoran's I",
#     # colors = brewer.pal(n = 11, name = "RdBu") %>% rev(),
#     colors = color_palette,
#     values = scales::rescale(c(-0.6, -0.3, 0, 0.3, 0.6)),
#     limits = c(-0.8, 0.8),
#     na.value = "grey"
#   ) +
#   theme(
#     panel.grid = element_blank(),
#     strip.placement = "outside",
#     strip.background.y = element_rect(fill = "grey95", color = "grey"),
#     strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
#     
#     axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
#     axis.title.x = element_blank(),
#     
#     axis.text.y = element_text(hjust = 1),
#     axis.title.y = element_blank(),
#     
#     legend.position = "right",
#     legend.box = "vertical",
#     legend.title = element_text(size = rel(0.8)),
#     legend.text = element_text(size = rel(0.8)),
#     text = element_text(size = 12),
#     
#     # 添加heatmap特定的美化选项
#     panel.spacing = unit(0.1, "lines"),  # 减少facet之间的间距
#     panel.border = element_rect(fill = NA, color = "grey80"),  # 添加边框区分facet
#   )
# p
# my_save_pdf(p, 'plots_manuscript/ravi_bimoran_tfa_mge.pdf', width = 8, height = 9)

