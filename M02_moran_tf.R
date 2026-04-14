setwd("~/Desktop/Pitt/Project_zhang_gbm_v2/")
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggpubr) 
library(ggrepel)
library(matrixStats)
library(patchwork)
library(RColorBrewer)
library(stringr)
library(tidyverse)
library(viridis)
source('aux_correlation.R')

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Part 1: RAVI data - all TFs
# ========== ========== ========== ========== ========== ========== ========== ==========

# Load Ravi data
moran_tfa_ravi = read.csv('dataframes/morani_ravi/tfa_per_sample.csv', row.names=1)
moran_mrna_ravi = read.csv('dataframes/morani_ravi/mrna_tf_per_sample.csv', row.names=1)

# Find common TFs between mRNA and TFA for Ravi data
common_tfs_ravi <- intersect(rownames(moran_mrna_ravi), rownames(moran_tfa_ravi))

# Filter to common TFs and replace NAs with 0
moran_tfa_ravi_fil = moran_tfa_ravi[common_tfs_ravi, ]
moran_mrna_ravi_fil = moran_mrna_ravi[common_tfs_ravi, ]
moran_tfa_ravi_fil[is.na(moran_tfa_ravi_fil)] <- 0
moran_mrna_ravi_fil[is.na(moran_mrna_ravi_fil)] <- 0

ht1 = Heatmap(moran_tfa_ravi_fil %>% as.matrix(),
              col = colorRamp2(seq(0, 1, length.out = 8),
                               brewer.pal(n=9, name = "Oranges")[1:8]),
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              column_names_side = "top",
              show_heatmap_legend = TRUE,
              show_row_names = FALSE,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              name = "TF\nActivity\nMoran's I")

ht2 = Heatmap(moran_mrna_ravi_fil %>% as.matrix(),
              col = colorRamp2(seq(0, 0.25, length.out = 8),
                               brewer.pal(n=9, name = "Blues")[1:8]),
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              column_names_side = "top",
              show_heatmap_legend = TRUE,
              show_row_names = FALSE,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              name = "mRNA\nExpression\nMoran's I")

p = ht1 + ht2
p

# my_save_pdf(p, 'plots_manuscript/ravi_moran_tf.pdf', width = 7.5, height = 7.5)

# Save the column order from Ravi data (TFs order from ht1)
ravi_tf_order <- rownames(moran_tfa_ravi_fil)[row_order(ht1)]

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Part 2: REN + SONPATKI data - all TFs
# ========== ========== ========== ========== ========== ========== ========== ==========

# Load Ren and Sonpatki data and merge
moran_tfa_ren <- merge(
  read.csv('dataframes/morani_ren/tfa_per_sample.csv', row.names=1),
  read.csv('dataframes/morani_sonpatki/tfa_per_sample.csv', row.names=1),
  by = "row.names", 
  all = TRUE
) %>% column_to_rownames("Row.names")

moran_mrna_ren <- merge(
  read.csv('dataframes/morani_ren/mrna_tf_per_sample.csv', row.names=1),
  read.csv('dataframes/morani_sonpatki/mrna_tf_per_sample.csv', row.names=1),
  by = "row.names", 
  all = TRUE
) %>% column_to_rownames("Row.names")


# Find common TFs between Ren+Sonpatki mRNA and TFA data, AND also present in Ravi order
common_tfs_ren <- intersect(
  intersect(rownames(moran_mrna_ren), rownames(moran_tfa_ren)),
  ravi_tf_order
)

# Filter to common TFs
moran_tfa_ren_fil = moran_tfa_ren[common_tfs_ren, ]
moran_mrna_ren_fil = moran_mrna_ren[common_tfs_ren, ]

# Replace NAs with 0
moran_tfa_ren_fil[is.na(moran_tfa_ren_fil)] <- 0
moran_mrna_ren_fil[is.na(moran_mrna_ren_fil)] <- 0

# Reorder columns (TFs) to match Ravi order
# First, find which TFs from ravi_tf_order are present in our filtered data
tfs_in_order <- ravi_tf_order[ravi_tf_order %in% common_tfs_ren]

# Reorder the dataframes
moran_tfa_ren_fil <- moran_tfa_ren_fil[tfs_in_order, ]
moran_mrna_ren_fil <- moran_mrna_ren_fil[tfs_in_order, ]

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Part 3: scatter plot
# ========== ========== ========== ========== ========== ========== ========== ==========

# row_medians_ravi <- rowMedians(as.matrix(moran_pwa_ravi_fil[, sapply(moran_pwa_ravi_fil, is.numeric)]),
#                                na.rm = TRUE)
# row_medians_ren <- rowMedians(as.matrix(moran_pwa_ren_fil[, sapply(moran_pwa_ren_fil, is.numeric)]),
#                               na.rm = TRUE)
# df_moran <- data.frame(
#   ravi = row_medians_ravi,
#   ren = row_medians_ren,
#   row_name = rownames(moran_pwa_ravi_fil)
# ) 

df_ravi = read.csv('dataframes/median_moran_ravi.csv', row.names = 1)
df_ren = read.csv('dataframes/median_moran_ren.csv', row.names = 1)

df_moran <- merge(
  df_ravi[, c("row_name", "TFA_m")], 
  df_ren[, c("row_name", "TFA_m")], 
  by = "row_name",
  suffixes = c("_ravi", "_ren")
) 
colnames(df_moran) <- c("row_name", "ravi", "ren")
df_moran_clean <- df_moran %>%
  filter(ravi != 0 & ren != 0)

df_moran_dotplot <- df_moran_clean %>%
  mutate(
    label = case_when(
      ravi > 0.98 & ren > 0.98 ~ row_name,
      TRUE ~ ""
    )
  )

s = cor.test(df_moran_clean$ravi, df_moran_clean$ren, method = "spearman")
print(s$estimate)

p <- ggplot(df_moran_dotplot, aes(x = ravi, y = ren)) +
  geom_abline(intercept = 0, slope = 1, color = "#888888", linetype = "dashed", size = 1) +
  annotate("text", x = 0.75, y = 0.0,
           label = expression(paste("Spearman ", rho, " = ", 0.9436)),
           size = 4, hjust = 0, vjust = 0.5) +
  geom_point(size = 2, color = "#313695", shape = 21, fill = "#4575B4") +
  geom_label_repel(
    aes(label = label),
    size = 4, max.overlaps = Inf, box.padding = 0.5,
    label.padding = 0.1, label.size = 0.2,
    segment.color = "#666666",
    segment.size = 0.5,
    min.segment.length = 0,
    force = 20, force_pull = 1,
  # direction = "both", nudge_x = 0.02, nudge_y = 0.02
  ) +
  xlim(0, 1.3) + ylim(0, 1.3) + 
  labs(x = "Discovery cohort (Median Moran's I)", 
       y = "Validation cohort (Median Moran's I)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    panel.grid = element_blank()
  )

p
my_save_pdf(p, 'plots_manuscript/dot_moran_tf.pdf', width = 6, height = 6)

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Part 4: RAVI data - filtered TFs
# ========== ========== ========== ========== ========== ========== ========== ==========

moran_tfa_ravi_subset = filter_dataframe_by_threshold(
  moran_tfa_ravi, thre_val = 0.9, thre_sam = 18)
common_tfs = rownames(moran_tfa_ravi_subset)
moran_tfa_ravi_subset[is.na(moran_tfa_ravi_subset)] <- 0
moran_mrna_ravi_subset = moran_mrna_ravi_fil[common_tfs, ]

# ht1_subset = Heatmap(moran_tfa_ravi_subset %>% as.matrix(),
#                      col = colorRamp2(seq(0.8, 1, length.out = 8),
#                                       brewer.pal(n=9, name = "Oranges")[1:8]),
#                      cluster_rows = TRUE,
#                      cluster_columns = FALSE,
#                      show_column_dend = FALSE,
#                      show_row_dend = FALSE,
#                      column_names_side = "top",
#                      row_names_side = "left",
#                      show_heatmap_legend = TRUE,
#                      show_row_names = TRUE,
#                      row_names_gp = gpar(fontsize = 8),  
#                      column_names_gp = gpar(fontsize = 8), 
#                      name = "TF\nActivity\nMoran's I")
# ht1_subset
# 
# ht2_subset = Heatmap(moran_mrna_ravi_subset %>% as.matrix(),
#                      col = colorRamp2(seq(0, 0.25, length.out = 8),
#                                       brewer.pal(n=9, name = "Blues")[1:8]),
#                      cluster_rows = TRUE,
#                      cluster_columns = FALSE,
#                      show_column_dend = FALSE,
#                      show_row_dend = FALSE,
#                      column_names_side = "top",
#                      row_names_side = "left",
#                      show_heatmap_legend = TRUE,
#                      show_row_names = TRUE,
#                      row_names_gp = gpar(fontsize = 8),  
#                      column_names_gp = gpar(fontsize = 8), 
#                      name = "mRNA\nExpression\nMoran's I")
# ht2_subset
# 
# p = ht1_subset + ht2_subset
# p
# 
# my_save_pdf(p, 'plots_manuscript/ravi_moran_tf_subset1.pdf', width = 7.5, height = 6)

# ========== ========== ========== ==========

tf_groups <- list(
  "Stress" = c("ATF5", "NFE2L2"),
  "Hypoxia\nMetabolic" = c("SREBF2", "HSF2", "EPAS1"),
  "Signaling" = c("LEF1", "SMAD4", "SNAI2"),
  "Immune\nInflammatory" = c("IKZF1", "TAL1", "MAFB", "MAFG", "STAT3", "STAT6", "NR3C1"),
  "Lineage\nDevelopmental" = c("PAX6", "SOX2", "SOX10", "PRDM1", "RFX3", "MEF2C", "LHX2", "RARA", "HES1"),
  "ETS Family" = c("ELF2", "ETV4", "ETV6"),
  "ZNF Family" = c("ZNF12", "ZNF217", "ZNF236", "ZNF250", "ZNF281", "ZNF711", "ZNF84"),
  "Other" = c("CREB3", "CUX1", "DBP", "HBP1", "MBD2", "MZF1", "NFYC", "TP53", "ZFHX3")
)
desired_order <- unlist(tf_groups)
tf_group_df <- data.frame(
  tf_gene = desired_order,
  cell_type = factor(rep(names(tf_groups), times = sapply(tf_groups, length)),
                     levels = names(tf_groups)))

tfa_long <- moran_tfa_ravi_subset[desired_order, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("tf_gene") %>%
  pivot_longer(-tf_gene, names_to = "sample", values_to = "moran_i") %>%
  left_join(tf_group_df, by = "tf_gene") %>%
  mutate(tf_gene = factor(tf_gene, levels = desired_order))
mrna_long <- moran_mrna_ravi_subset[desired_order, , drop = FALSE] %>%
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
    name = "TF\nActivity\nMoran's I",
    colors = brewer.pal(n = 9, name = "Oranges")[1:8],
    values = scales::rescale(seq(0.8, 1, length.out = 8)),
    limits = c(0.8, 1), 
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
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    text = element_text(size = 10),
    panel.spacing = unit(0.1, "lines"),
    panel.border = element_rect(fill = NA, color = "grey80")
  )

p_tfa

# my_save_pdf(p_tfa, 'plots_manuscript/ravi_moran_tf_heatmap1.pdf', width = 5.4, height = 6)

# Create right plot (mRNA Expression) - remove y-axis AND strip labels
p_mrna <- ggplot(mrna_long, aes(x = sample, y = tf_gene)) +
  geom_tile(aes(fill = moran_i)) +
  facet_grid(cell_type ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_classic() +
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(
    name = "mRNA\nExpression\nMoran's I",
    colors = brewer.pal(n = 9, name = "Blues")[1:8],
    values = scales::rescale(seq(0, 0.25, length.out = 8)),
    limits = c(0, 0.20),
    na.value = "white"
  ) +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background.y = element_blank(),  # Remove strip background
    strip.text.y.left = element_blank(),    # Remove strip text
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    text = element_text(size = 10),
    panel.spacing = unit(0.1, "lines"),
    panel.border = element_rect(fill = NA, color = "grey80")
  )

# Combine the plots
combined_plot <- p_tfa + p_mrna +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(combined_plot)
my_save_pdf(combined_plot, 'plots_manuscript/ravi_moran_tf_heatmap1.pdf', width = 9, height = 6)

# ========== ========== ========== ==========

moran_mrna_ravi_subset = filter_dataframe_by_threshold(
  moran_mrna_ravi, thre_val = 0.05, thre_sam = 18)
common_mrnas = setdiff(rownames(moran_mrna_ravi_subset), common_tfs)
moran_mrna_ravi_subset[is.na(moran_mrna_ravi_subset)] <- 0
moran_tfa_ravi_subset = moran_tfa_ravi_fil[common_mrnas, ]

ht1_subset = Heatmap(moran_tfa_ravi_subset %>% as.matrix(),
                     col = colorRamp2(seq(0, 1, length.out = 8),
                                      brewer.pal(n=9, name = "Oranges")[1:8]),
                     cluster_rows = TRUE,
                     cluster_columns = FALSE,
                     show_column_dend = FALSE,
                     show_row_dend = FALSE,
                     column_names_side = "top",
                     row_names_side = "left",
                     show_heatmap_legend = TRUE,
                     show_row_names = TRUE,
                     row_names_gp = gpar(fontsize = 8),  
                     column_names_gp = gpar(fontsize = 8), 
                     name = "TF\nActivity\nMoran's I")
ht1_subset

ht2_subset = Heatmap(moran_mrna_ravi_subset[common_mrnas, ] %>% as.matrix(),
                     col = colorRamp2(seq(0, 0.25, length.out = 8),
                                      brewer.pal(n=9, name = "Blues")[1:8]),
                     cluster_rows = TRUE,
                     cluster_columns = FALSE,
                     show_column_dend = FALSE,
                     show_row_dend = FALSE,
                     column_names_side = "top",
                     row_names_side = "left",
                     show_heatmap_legend = TRUE,
                     show_row_names = TRUE,
                     row_names_gp = gpar(fontsize = 8),  
                     column_names_gp = gpar(fontsize = 8), 
                     name = "mRNA\nExpression\nMoran's I")
ht2_subset

p = ht1_subset + ht2_subset
p

my_save_pdf(p, 'plots_manuscript/ravi_moran_tf_subset2.pdf', width = 12.5, height = 2.5)

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Part 5: REN + SONPATKI data - filtered TFs
# ========== ========== ========== ========== ========== ========== ========== ==========

moran_tfa_ren_subset = filter_dataframe_by_threshold(
  moran_tfa_ren, thre_val = 0.8, thre_sam = 13)
common_tfs_ren = rownames(moran_tfa_ren_subset) %>% intersect(rownames(moran_tfa_ravi_subset))
moran_tfa_ren_subset = moran_tfa_ren[common_tfs_ren, ]
moran_tfa_ren_subset[is.na(moran_tfa_ren_subset)] <- 0
moran_mrna_ren_subset = moran_mrna_ren_fil[common_tfs_ren, ]

tf_groups_ren = list(
  "Stress" = c("ATF5", "NFE2L2"),
  "Hypoxia\nMetabolic" = c("SREBF2", "HSF2", "EPAS1"),
  "Signaling" = c("LEF1", "SMAD4", "SNAI2"),
  "Immune\nInflammatory" = c("IKZF1", "MAFB", "MAFG", "STAT3", "STAT6", "NR3C1"),  # 删除了 TAL1
  "Lineage\nDevelopmental" = c("PAX6", "SOX2", "SOX10", "PRDM1", "RFX3", "MEF2C", "LHX2", "RARA"),  # 删除了 HES1
  "ETS Family" = c("ELF2", "ETV4"),  # 删除了 ETV6
  "ZNF Family" = c("ZNF12", "ZNF217", "ZNF236", "ZNF250", "ZNF281", "ZNF711", "ZNF84"),
  "Other" = c("CREB3", "CUX1", "DBP", "HBP1", "MBD2", "MZF1", "NFYC", "TP53", "ZFHX3")
)
desired_order_ren <- unlist(tf_groups_ren)
tf_group_df_ren <- data.frame(
  tf_gene = desired_order_ren,
  cell_type = factor(rep(names(tf_groups_ren), times = sapply(tf_groups_ren, length)),
                     levels = names(tf_groups_ren)))

tfa_long_ren <- moran_tfa_ren_subset[desired_order_ren, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("tf_gene") %>%
  pivot_longer(-tf_gene, names_to = "sample", values_to = "moran_i") %>%
  left_join(tf_group_df_ren, by = "tf_gene") %>%
  mutate(tf_gene = factor(tf_gene, levels = desired_order_ren))
mrna_long_ren <- moran_mrna_ren_subset[desired_order_ren, , drop = FALSE] %>%
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
    name = "TF\nActivity\nMoran's I",
    colors = brewer.pal(n = 9, name = "Oranges")[1:8],
    values = scales::rescale(seq(0.8, 1, length.out = 8)),
    limits = c(0.8, 1), 
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
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    text = element_text(size = 10),
    panel.spacing = unit(0.1, "lines"),
    panel.border = element_rect(fill = NA, color = "grey80")
  )

p_tfa_ren

# my_save_pdf(p_tfa, 'plots_manuscript/ravi_moran_tf_heatmap1.pdf', width = 5.4, height = 6)

# Create right plot (mRNA Expression) - remove y-axis AND strip labels
p_mrna_ren <- ggplot(mrna_long_ren, aes(x = sample, y = tf_gene)) +
  geom_tile(aes(fill = moran_i)) +
  facet_grid(cell_type ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_classic() +
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(
    name = "mRNA\nExpression\nMoran's I",
    colors = brewer.pal(n = 9, name = "Blues")[1:8],
    values = scales::rescale(seq(0, 0.25, length.out = 8)),
    limits = c(0, 0.20),
    na.value = "white"
  ) +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background.y = element_blank(),  # Remove strip background
    strip.text.y.left = element_blank(),    # Remove strip text
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    text = element_text(size = 10),
    panel.spacing = unit(0.1, "lines"),
    panel.border = element_rect(fill = NA, color = "grey80")
  )

# Combine the plots
combined_plot <- p_tfa_ren + p_mrna_ren +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(combined_plot)
my_save_pdf(combined_plot, 'plots_manuscript/ren_moran_tf_heatmap1.pdf', width = 7.5, height = 6)

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Part 6: Human protein atlas
# ========== ========== ========== ========== ========== ========== ========== ==========

library(BiocStyle)
library(clusterProfiler)
library(HPAanalyze)
library(hpar)      
library(httr)    
library(org.Hs.eg.db)  
library(patchwork)
library(tibble)   
library(tidyr)
library(xml2)

# https://bioconductor.posit.co/packages/devel/bioc/vignettes/HPAanalyze/inst/doc/b_HPAanalyze_indepth.html

fname = 'plots_gbm_hpa_images/glioma_pathology.rds'
if (file.exists(fname)) {
  glioma_pat_data <- readRDS(fname)
} else {
  downloadedData <- hpaDownload(downloadList=c("pathology"))
  glioma_pat_data <- hpaSubset(data=downloadedData, targetCancer='glioma')
  saveRDS(glioma_pat_data, file = fname)
}

glioma_pat_data$pathology_renamed <- glioma_pat_data$pathology %>%
  rename_with(~ c("High", "Medium", "Low", "Not Detected"), 
              .cols = c("high", "medium", "low", "not_detected"))

fname = 'plots_gbm_hpa_images/glioma_normal.rds'
if (file.exists(fname)) {
  glioma_normal_data <- readRDS(fname)
} else {
  downloadedData <- hpaDownload(downloadList=c("normal_tissue"))
  glioma_normal_data <- hpaSubset(data=downloadedData, 
                                  targetTissue = "Cerebral cortex")
  saveRDS(glioma_normal_data, file = fname)
}

glioma_normal_data$normal_renamed <- glioma_normal_data$normal_tissue %>%
  mutate(level = case_when(
    level == "high" ~ "High",
    level == "medium" ~ "Medium",
    level == "low" ~ "Low",
    level == "Not detected" ~ "Not Detected", 
    TRUE ~ level
  ))

# ========== ========== ==========

plot_protein_atlas <- function(tf_list){
  glioma_pat_subset <- glioma_pat_data$pathology_renamed %>%
    filter(gene %in% tf_list) %>%  # 仅筛选tf_list中的基因，无其他过滤
    pivot_longer(
      cols = c(High, Medium, Low, `Not Detected`),
      names_to = "expression_level",
      values_to = "count"
    ) %>%
    group_by(gene) %>%  
    mutate( percent = count / sum(count) * 100 ) %>%
    ungroup()  %>%
    mutate(expression_level = factor(expression_level,
                                     levels = c("Not Detected", "Low", "Medium", "High")))
  
  p = ggplot(glioma_pat_subset, 
             aes(x = gene, y = percent, fill = expression_level)) +
    geom_bar(position = "fill", stat = "identity", width = 0.8) + 
    scale_fill_manual(values = brewer.pal(n=4, name = "YlOrBr"),
                      name = 'Level') + 
    theme_bw() + 
    labs( x = '',
          y = 'Patient Proportion') + 
    theme(
      # 保留之前的隐藏设置
      axis.line = element_blank(),        # 隐藏坐标轴轴线
      axis.ticks = element_blank(),       # 隐藏刻度线
      panel.grid = element_blank(),       # 隐藏背景虚线
      panel.border = element_blank(),     # 隐藏外边框（可选，根据需求保留/删除）
      
      # 核心修改：调整x轴标签位置，贴近柱状图
      axis.text.x = element_text(
        angle = 90,        # 标签垂直显示
        hjust = 0.5,       # 水平居中
        vjust = 0.5,       # 垂直偏移（减小值让标签更贴近bar，默认0.5）
        margin = margin(t = -5, unit = "pt")  # 向上偏移5个点（t=top，负数=向上）
      ),
      text = element_text(size=15)
    )
  
  glioma_normal_subset <- glioma_normal_data$normal_renamed %>%
    filter(gene %in% common_tfs) %>%
    group_by(gene) %>%  
    ungroup() %>%
    mutate(level = factor(level,
                          levels = c("Not Detected", "Low", "Medium", "High")))
  
  q = ggplot(glioma_normal_subset, 
             aes(x = `cell_type`, y = gene, fill = level)) +
    geom_tile(color = "white", size = 0.2) +
    scale_fill_manual(values = brewer.pal(n=4, name = "YlOrBr"),
                      name = 'Level',
                      guide = guide_legend(ncol = 4)) +  # 设置legend为4列
    theme_bw() + coord_flip() +
    labs( y = '',
          x = 'Cell Type') + 
    theme(
      # 保留之前的隐藏设置
      axis.line = element_blank(),        # 隐藏坐标轴轴线
      axis.ticks = element_blank(),       # 隐藏刻度线
      panel.grid = element_blank(),       # 隐藏背景虚线
      panel.border = element_blank(),     # 隐藏外边框（可选，根据需求保留/删除）
      
      # 核心修改：调整x轴标签位置，贴近柱状图
      axis.text.x = element_text(
        angle = 90,        # 标签垂直显示
        hjust = 0.5,       # 水平居中
        vjust = 0.5,       # 垂直偏移（减小值让标签更贴近bar，默认0.5）
        margin = margin(t = 0, unit = "pt")  # 向上偏移5个点（t=top，负数=向上）
      ),
      text = element_text(size=15),
      
      # 调整legend位置到图2下面
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.margin = margin(t = -10, unit = "pt"),
      legend.spacing.x = unit(0.5, "cm")  # 设置legend项之间的间距
    )
  
  # 移除图1的legend（如果图1也有legend）
  p_no_legend <- p + theme(legend.position = "none")
  
  # 组合图形
  combined_plot <- p_no_legend / q + 
    plot_layout(heights = c(20/3, 10/3))
  
  return(combined_plot)
}

# tf_list = union(common_tfs, common_mrnas)
combined_plot = plot_protein_atlas(common_tfs)
combined_plot
my_save_pdf(combined_plot, 'plots_manuscript/ravi_hpa_tf.pdf', width = 9, height = 5.4)


