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
#           Part 1: RAVI data - generate column order
# ========== ========== ========== ========== ========== ========== ========== ==========

# Load Ravi pathway data
moran_pwa_ravi = read.csv('dataframes/morani_ravi/pwa_per_sample.csv', row.names=1)

# Process Ravi data
moran_pwa_ravi_fil = moran_pwa_ravi
rownames(moran_pwa_ravi_fil) = hallmarks[rownames(moran_pwa_ravi)]
moran_pwa_ravi_fil[is.na(moran_pwa_ravi_fil)] <- 0

# Create heatmap for Ravi data (with clustering to determine order)
ht3 = Heatmap(moran_pwa_ravi_fil %>% as.matrix(),
              col = colorRamp2(seq(0, 1, length.out = 8),
                               brewer.pal(n=9, name = "Greens")[1:8]),
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              column_names_side = "top",
              show_heatmap_legend = TRUE,
              show_row_names = TRUE,
              row_names_gp = gpar(fontsize = 8),  
              column_names_gp = gpar(fontsize = 8), 
              name = "Pathway\nActivity\nMoran's I")

ht3

# Save Ravi plot
my_save_pdf(ht3, 'plots_manuscript/ravi_moran_pw.pdf', width = 6, height = 7.5)

# Save the column order from Ravi data (pathway order)
ravi_pw_order <- rownames(moran_pwa_ravi_fil)[row_order(ht3)]

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Part 2: REN + SONPATKI data - use same column order
# ========== ========== ========== ========== ========== ========== ========== ==========

# Load Ren and Sonpatki pathway data and merge
moran_pwa_ren <- merge(
  read.csv('dataframes/morani_ren/pwa_per_sample.csv', row.names=1),
  read.csv('dataframes/morani_sonpatki/pwa_per_sample.csv', row.names=1),
  by = "row.names", 
  all = TRUE
) %>% column_to_rownames("Row.names")

# Process Ren+Sonpatki data
moran_pwa_ren_fil = moran_pwa_ren
rownames(moran_pwa_ren_fil) = hallmarks[rownames(moran_pwa_ren)]
moran_pwa_ren_fil[is.na(moran_pwa_ren_fil)] <- 0

# Reorder the dataframe columns to match Ravi order
moran_pwa_ren_fil <- moran_pwa_ren_fil[ravi_pw_order, ]

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Part 3: scatter plot
# ========== ========== ========== ========== ========== ========== ========== ==========

row_medians_ravi <- rowMedians(as.matrix(moran_pwa_ravi_fil[, sapply(moran_pwa_ravi_fil, is.numeric)]),
                               na.rm = TRUE)
row_medians_ren <- rowMedians(as.matrix(moran_pwa_ren_fil[, sapply(moran_pwa_ren_fil, is.numeric)]),
                              na.rm = TRUE)
df_moran <- data.frame(
  ravi = row_medians_ravi,
  ren = row_medians_ren,
  row_name = rownames(moran_pwa_ravi_fil)
) 

df_moran_dotplot <- df_moran %>%
  mutate(
    label = case_when(
      ravi > 0.95 & ren > 0.95 ~ row_name,
      TRUE ~ ""
    )
  )

s = cor.test(df_moran$ravi, df_moran$ren, method = "spearman")
print(s$estimate)

p <- ggplot(df_moran_dotplot, aes(x = ravi, y = ren)) +
  geom_abline(intercept = 0, slope = 1, color = "#888888", linetype = "dashed", size = 1) +
  annotate("text", x = 0.95, y = 0.85,
           label = expression(paste("Spearman ", rho, " = ", -0.0159)),
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
  xlim(0.84, 1.02) + ylim(0.84, 1.02) + 
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
my_save_pdf(p, 'plots_manuscript/dot_moran_pw.pdf', width = 6, height = 6)

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Part 4: RAVI data - filtered PWs
# ========== ========== ========== ========== ========== ========== ========== ==========

moran_pwa_ravi_subset = filter_dataframe_by_threshold(
  moran_pwa_ravi_fil, thre_val = 0.9, thre_sam = 18)
common_pws = rownames(moran_pwa_ravi_subset)
moran_pwa_ravi_subset[is.na(moran_pwa_ravi_subset)] <- 0

# Create heatmap for Ravi data (with clustering to determine order)
ht3 = Heatmap(moran_pwa_ravi_subset %>% as.matrix(),
              col = colorRamp2(seq(0.8, 1, length.out = 8),
                               brewer.pal(n=9, name = "Greens")[1:8]),
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              column_names_side = "top",
              row_names_side = "left",
              show_heatmap_legend = TRUE,
              show_row_names = TRUE,
              row_names_gp = gpar(fontsize = 10),  
              column_names_gp = gpar(fontsize = 10), 
              name = "Pathway\nActivity\nMoran's I")

ht3

# Save Ravi plot
my_save_pdf(ht3, 'plots_manuscript/ravi_moran_pw_subset.pdf', width = 7, height = 6)

