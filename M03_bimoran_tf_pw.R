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
    parts <- str_split(x, "_", n = 3)[[1]] 
    if(length(parts) >= 2 && paste(parts[2], parts[3], sep = "_") %in% names(hallmarks)) {
      pathway_name <- paste(parts[2], parts[3], sep = "_")
      readable_name <- hallmarks[[pathway_name]]
      return(paste(parts[1], readable_name, sep = " - "))
    } else if(parts[2] %in% names(hallmarks)) {
      readable_name <- hallmarks[[parts[2]]]
      return(paste(parts[1], readable_name, sep = " - "))
    } else {
      return(x)
    }
  })
  rownames(df) <- new_rownames
  return(df)
}

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Bivariate: TF x PW
# ========== ========== ========== ========== ========== ========== ========== ==========

bivariate_tfa_pwa_ravi = read.csv('dataframes/morani_ravi/bivariate_tfa_pwa_per_sample.csv', row.names=1) %>% 
  rename_rownames()

bivariate_tfa_pwa_ren = merge(
  read.csv('dataframes/morani_ren/bivariate_tfa_pwa_per_sample.csv', row.names=1),
  read.csv('dataframes/morani_sonpatki/bivariate_tfa_pwa_per_sample.csv', row.names=1),
  by = "row.names", 
  all = TRUE
) %>% column_to_rownames("Row.names") %>% 
  rename_rownames()

# Find common TF-pathway pairs
common_tfs <- intersect(rownames(bivariate_tfa_pwa_ravi), rownames(bivariate_tfa_pwa_ren))
bivariate_tfa_pwa_ravi_fil = bivariate_tfa_pwa_ravi[common_tfs, ]
bivariate_tfa_pwa_ren_fil = bivariate_tfa_pwa_ren[common_tfs, ]

# Replace NAs with 0
bivariate_tfa_pwa_ravi_fil[is.na(bivariate_tfa_pwa_ravi_fil)] <- 0
bivariate_tfa_pwa_ren_fil[is.na(bivariate_tfa_pwa_ren_fil)] <- 0

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Dotplot / Density plot
# ========== ========== ========== ========== ========== ========== ========== ==========

row_medians_ravi <- rowMedians(
  as.matrix(bivariate_tfa_pwa_ravi_fil[, sapply(bivariate_tfa_pwa_ravi_fil, is.numeric)]), na.rm = TRUE)
row_medians_ren <- rowMedians(
  as.matrix(bivariate_tfa_pwa_ren_fil[, sapply(bivariate_tfa_pwa_ren_fil, is.numeric)]), na.rm = TRUE)

df_moran <- data.frame(
  ravi = row_medians_ravi,
  ren = row_medians_ren,
  row_name = rownames(bivariate_tfa_pwa_ravi_fil)
) 

df_moran_clean <- df_moran %>%
  filter(ravi != 0 & ren != 0)
s = cor.test(df_moran_clean$ravi, df_moran_clean$ren, method = "spearman")
print(s$estimate)

df_moran_dotplot <- df_moran_clean %>%
  mutate(
    label = case_when(
      ravi > 0.6 & ren > 0.6 ~ row_name,
      ravi < -0.6 & ren < -0.6 ~ row_name,
      TRUE ~ ""
    )
  )

p <- ggplot(df_moran_dotplot, aes(x = ravi, y = ren)) +
  geom_abline(intercept = 0, slope = 1, color = "#888888", linetype = "dashed", size = 1) +
  annotate("text", x = 0.5, y = -0.8,
           label = expression(paste("Spearman ", rho, " = ", 0.5508)),
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
  xlim(-0.8, 1.3) + ylim(-0.8, 1.3) + 
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

my_save_pdf(p, 'plots_manuscript/dot_bimoran_tf_pw.pdf', width = 8, height = 8)

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Heatmap
# ========== ========== ========== ========== ========== ========== ========== ==========

bivariate_tfa_pwa_ravi_subset = filter_dataframe_by_threshold(
  bivariate_tfa_pwa_ravi_fil, thre_val = 0.4, thre_sam = 18)
common_tfs = rownames(bivariate_tfa_pwa_ravi_subset)
common_tfs

tf_groups <- list(
  "ATF5" = c("ATF5 - Apical Surface", "ATF5 - Heme Metabolism", 
             "ATF5 - KRAS Signaling Down", "ATF5 - Notch Signaling", 
             "ATF5 - Peroxisome"),
  "FLI1" = c("FLI1 - Cholesterol Homeostasis", "FLI1 - DNA Repair"),
  "HIF1A" = c("HIF1A - Hypoxia", "HIF1A - mTORC1 Signaling"),
  "LHX2" = c("LHX2 - TNFA Signaling via NFKB"),
  "MAZ" = c("MAZ - Androgen Response", "MAZ - Apical Junction", 
            "MAZ - Apical Surface", "MAZ - Coagulation", 
            "MAZ - Complement", "MAZ - DNA Repair", 
            "MAZ - E2F Targets", "MAZ - ER Early", 
            "MAZ - ER Late", "MAZ - Glycolysis", 
            "MAZ - Hedgehog Signalling", "MAZ - Heme Metabolism", 
            "MAZ - IL2-STAT5 Signaling", "MAZ - KRAS Signaling Down", 
            "MAZ - KRAS Signaling Up", "MAZ - Mitotic Spindle", 
            "MAZ - MYC Targets v1", "MAZ - Notch Signaling", 
            "MAZ - Oxidative phosphorylation", "MAZ - Peroxisome", 
            "MAZ - PI3K/AKT/mTOR Signaling", "MAZ - Protein Secretion", 
            "MAZ - TGFB Signaling", "MAZ - UV Response Down", 
            "MAZ - UV Response Up", "MAZ - Xenobiotic Metabolism"),
  "NR3C1" = c("NR3C1 - Hypoxia"),
  "RFX3" = c("RFX3 - Allograft Rejection"),
  "SMAD3" = c("SMAD3 - Hypoxia"),
  "SNAI2" = c("SNAI2 - EMT", "SNAI2 - IL6/JAK/STAT3 Signaling"),
  "SOX2" = c("SOX2 - Angiogenesis", "SOX2 - EMT"),
  "STAT2" = c("STAT2 - INFA Response"),
  "ZNF281" = c("ZNF281 - Apical Surface", "ZNF281 - UV Response Up")
)

desired_order <- unique(unlist(tf_groups))
tf_group_df <- data.frame(
  tf_gene = desired_order,
  cell_type = factor(rep(names(tf_groups), times = sapply(tf_groups, length)),
                     levels = names(tf_groups))
)

tfa_long <- bivariate_tfa_pwa_ravi_subset[desired_order, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("tf_gene") %>%
  pivot_longer(-tf_gene, names_to = "sample", values_to = "moran_i") %>%
  left_join(tf_group_df, by = "tf_gene") %>%
  mutate(
    tf_gene_clean = sub("^[A-Z0-9]+ - ", "", tf_gene)
  ) %>%
  mutate(
    tf_gene = factor(tf_gene_clean, levels = unique(tf_gene_clean))
  )

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

my_save_pdf(p_tfa, 'plots_manuscript/ravi_bimoran_tfa_pwa.pdf', width = 5, height = 7.5)

# ========== ========== ========== ========== ========== ========== ========== ==========
#           Validation
# ========== ========== ========== ========== ========== ========== ========== ==========

bivariate_tfa_pwa_ren_subset = filter_dataframe_by_threshold(
  bivariate_tfa_pwa_ren_fil, thre_val = 0.3, thre_sam = 13)
common_tfs_ren = rownames(bivariate_tfa_pwa_ren_subset) %>% intersect(rownames(bivariate_tfa_pwa_ravi_subset))
common_tfs_ren

tf_groups_ren <- list(
  "ATF5" = c("ATF5 - Heme Metabolism", "ATF5 - Notch Signaling", "ATF5 - Peroxisome"),
  "FLI1" = c("FLI1 - DNA Repair"),
  "HIF1A" = c("HIF1A - Hypoxia"),
  "LHX2" = c("LHX2 - TNFA Signaling via NFKB"),
  "MAZ" = c("MAZ - Androgen Response", "MAZ - Apical Junction", "MAZ - Apical Surface", 
            "MAZ - DNA Repair", "MAZ - ER Early", "MAZ - ER Late", 
            "MAZ - Glycolysis", "MAZ - Hedgehog Signalling", "MAZ - Heme Metabolism", 
            "MAZ - IL2-STAT5 Signaling", "MAZ - KRAS Signaling Down", "MAZ - KRAS Signaling Up", 
            "MAZ - Mitotic Spindle", "MAZ - Notch Signaling", "MAZ - Oxidative phosphorylation", 
            "MAZ - Peroxisome", "MAZ - PI3K/AKT/mTOR Signaling", "MAZ - Protein Secretion", 
            "MAZ - UV Response Down", "MAZ - UV Response Up", "MAZ - Xenobiotic Metabolism"),
  "RFX3" = c("RFX3 - Allograft Rejection"),
  "SMAD3" = c("SMAD3 - Hypoxia"),
  "SNAI2" = c("SNAI2 - EMT", "SNAI2 - IL6/JAK/STAT3 Signaling"),
  "SOX2" = c("SOX2 - Angiogenesis"),
  "STAT2" = c("STAT2 - INFA Response")
)

desired_order_ren <- unlist(tf_groups_ren)
tf_group_df_ren <- data.frame(
  tf_gene = desired_order_ren,
  cell_type = factor(rep(names(tf_groups_ren), times = sapply(tf_groups_ren, length)),
                     levels = names(tf_groups_ren)))

tfa_long_ren <- bivariate_tfa_pwa_ren_subset[desired_order_ren, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("tf_gene") %>%
  pivot_longer(-tf_gene, names_to = "sample", values_to = "moran_i") %>%
  left_join(tf_group_df_ren, by = "tf_gene") %>%
  mutate(
    tf_gene_clean = sub("^[A-Z0-9]+ - ", "", tf_gene)
  ) %>%
  mutate(
    tf_gene = factor(tf_gene_clean, levels = unique(tf_gene_clean))
  )

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

my_save_pdf(p_tfa_ren, 'plots_manuscript/ren_bimoran_tfa_pwa.pdf', width = 5, height = 7)
