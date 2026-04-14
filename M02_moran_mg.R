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

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Moran's I dotplots
# ========== ========== ========== ========== ========== ========== ========== ==========

mrna_row_medians <- rowMedians(as.matrix(moran_mrna_ravi_fil[, sapply(moran_mrna_ravi_fil, is.numeric)]),
                               na.rm = TRUE)
tfa_row_medians <- rowMedians(as.matrix(moran_tfa_ravi_fil[, sapply(moran_tfa_ravi_fil, is.numeric)]),
                              na.rm = TRUE)
df_moran_dotplot <- data.frame(
  mRNA_m = mrna_row_medians,
  TFA_m = tfa_row_medians,
  row_name = names(mrna_row_medians)
) %>%
  mutate(
    label = case_when(
      mRNA_m > 0.1 ~ row_name,
      TFA_m > 0.98 ~ row_name,
      TRUE ~ ""
    )
  )
data.frame(
  mRNA_m = mrna_row_medians,
  TFA_m = tfa_row_medians,
  row_name = names(mrna_row_medians)
) %>% write.csv(file = 'dataframes/median_moran_ravi.csv')

p <- ggplot(df_moran_dotplot, aes(x = mRNA_m, y = TFA_m)) +
  geom_abline(intercept = 0, slope = 1, color = "#F46D43", linetype = "dashed", size = 1) +
  annotate("text", x = 0.35, y = 0.3,
           label = expression(italic(y) == italic(x)),
           size = 4, color = "#F46D43", hjust = 0, vjust = 0.5) +
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
  coord_cartesian(xlim = c(-0.05, 0.4), ylim = c(0, 1.2)) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.5", "0.75", "1")) +
  labs(x = "mRNA Moran's I", y = "TF Moran's I") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    panel.grid = element_blank()
  )

print(p)
my_save_pdf(p, 'plots_manuscript/dot_moran_marker.pdf', width = 6, height = 6.75)

