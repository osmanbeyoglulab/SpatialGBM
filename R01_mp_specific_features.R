setwd("~/Desktop/Pitt/Project_zhang_gbm_v2/")
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(RColorBrewer)
library(scales)
source('aux_correlation.R')
dir.create('plots_website/mp_specific')

FONTSIZE = 15
RELSIZE = 0.8

df_mp = read.csv('dataframes/metaprogram_proportion.csv', row.names=1)
# set1 = c("AC", "Chromatin.Reg", "MES", "MES.Ast", "MES.Hyp", "NPC", "OPC" , "Prolif.Metab")
# set2 = df_mp %>% rownames() %>% setdiff(set1)

# ========== ========== ========== ========== ==========

df_mp$MP <- rownames(df_mp)
df_long <- pivot_longer(df_mp, cols = -MP, names_to = "sample", values_to = "value")

p = ggplot(df_long, aes(x = sample, y = value, fill = MP)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = metaprogram_palette) +
  scale_y_continuous(expand = c(0, 0)) +  # This removes padding at y-axis limits
  labs(x = "Sample", y = "Proportion", fill = "Metaprogram") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    text = element_text(size = FONTSIZE)
  )
p

p = ggplot(df_long, aes(y = sample, x = value, fill = MP)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = metaprogram_palette) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Sample", x = "Proportion", fill = "Metaprogram") +
  theme_bw() +
  theme(
    legend.position = "top",           # Move legend to top
    legend.direction = "horizontal",   # Arrange items horizontally
    legend.key.size = unit(0.4, "cm"), # Adjust legend key size if needed
    legend.text = element_text(size = rel(0.8)), # Adjust text size if needed
    legend.title = element_blank(),
    axis.title.x = element_text(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    text = element_text(size = FONTSIZE)
  ) +
  guides(fill = guide_legend(ncol = 3))  # Set 4 columns in the legend

p

pdf('plots_manuscript/mp_barplot.pdf', width=7, height=8)
print(p)
dev.off()

# ========== ========== ========== ========== ==========

cluster_tfs <- function(data) {
  mat <- data %>%
    select(TF, sample, diff) %>%
    pivot_wider(names_from = sample, values_from = diff, values_fill = 0) %>%
    column_to_rownames("TF")
  if (nrow(mat) <= 1) { return(rownames(mat)) }
  hc <- hclust(dist(mat))
  return(rownames(mat)[hc$order])
}

factor_tfs <- function(data){
  tf_order_list <- data %>%
    group_by(subtype) %>% group_map(~ cluster_tfs(.x))
  
  subtype_names <- data %>%
    group_by(subtype) %>% group_keys() %>% pull(subtype)
  
  tf_subtype_levels <- map2(tf_order_list, subtype_names, ~ paste(.x, .y, sep = "|")) %>%
    unlist()
  
  data <- data %>%
    mutate(TF_subtype = paste(TF, subtype, sep = "|"),
           TF_subtype = factor(TF_subtype, levels = tf_subtype_levels))
  return(data)
}

# ========== ========== ========== ========== ==========

final_result_tfs <- read.csv("dataframes/mp_specific_tfs.csv", row.names = 1) %>% factor_tfs()
p <- ggplot(final_result_tfs, aes(x = TF_subtype, y = sample)) +
  geom_point(aes(size = neg_log_p_adj, color = diff)) +
  scale_color_gradientn(
    'Difference\nin Activity',
    colors = rev(brewer.pal(11, 'PiYG')),
    breaks = seq(-1.5, 1.5, 0.5),
    limits = c(-1.75, 1.75),
    # breaks = seq(-0.5, 0.5, 0.5), 
    # limits = c(-0.75, 0.75),
    oob = scales::oob_squish
  ) +
  scale_size_continuous(
    range = c(1, 5),
    name = "-log10(adj_p)"
  ) +
  scale_x_discrete(
    labels = function(x) sub("\\|.*", "", x)  # Show only TF names on x-axis
  ) +
  facet_grid(. ~ subtype, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(
    x = "Transcription Factors",
    y = "Samples",
  ) + 
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(colour = "grey", fill = "white"),
    strip.text.x = element_text(angle = 90),
    strip.text.y.left = element_text(face = 'bold', angle = 0),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = rel(RELSIZE)),
    legend.text = element_text(size = rel(RELSIZE)),
    legend.key.height = unit(1, "cm"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    text = element_text(size = FONTSIZE)
  )

print(p)
my_save_plot(p, 'plots_manuscript/mp_specific_tfs.pdf', width = 15, height = 8)

# ========== ========== ========== ========== ==========

final_result_pws <- read.csv("dataframes/mp_specific_pws.csv", row.names = 1) %>% rename(TF = 2) %>% factor_tfs()
p <- ggplot(final_result_pws, aes(x = TF_subtype, y = sample)) +
  geom_point(aes(size = neg_log_p_adj, color = diff)) +
  scale_color_gradientn(
    'Difference\nin Activity',
    colors = rev(brewer.pal(11, 'PuOr')),
    breaks = seq(-1.5, 1.5, 0.5), 
    limits = c(-1.75, 1.75),
    oob = scales::oob_squish
  ) +
  scale_size_continuous(
    range = c(1, 5),
    name = "-log10(adj_p)"
  ) +
  scale_x_discrete(
    labels = function(x) sub("\\|.*", "", x)  # Show only TF names on x-axis
  ) +
  facet_grid(. ~ subtype, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(x = "Pathways", y = "Samples") + 
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(colour = "grey", fill = "white"),
    strip.text.x = element_text(angle = 90),
    strip.text.y.left = element_text(face = 'bold', angle = 0),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = rel(RELSIZE)),
    legend.text = element_text(size = rel(RELSIZE)),
    legend.key.height = unit(1, "cm"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    text = element_text(size = FONTSIZE)
  )

print(p)
my_save_plot(p, 'plots_manuscript/mp_specific_pws.pdf', width = 15, height = 10)

# ========== ========== ========== ========== ==========

final_result_dcs <- read.csv("dataframes/mp_specific_dcs.csv", row.names = 1) %>% rename(TF = 2) %>% factor_tfs()
p <- ggplot(final_result_dcs, aes(x = TF_subtype, y = sample)) +
  geom_point(aes(size = neg_log_p_adj, color = diff)) +
  scale_color_gradientn(
    'Difference\nin Score',
    colors = brewer.pal(9, 'Oranges'),
    breaks = seq(0, 3, 1), 
    limits = c(0, 3.25),
    oob = scales::oob_squish
  ) +
  scale_size_continuous(
    range = c(1, 5),
    name = "-log10(adj_p)"
  ) +
  # scale_x_discrete(
  #   labels = function(x) sub("\\|[^|]*$", "", x)  # Keep everything before the last pipe
  # ) +
  scale_x_discrete(
    labels = function(x) sub("^[^|]*\\|([^|]*)\\|.*", "\\1", x)
  ) +
  facet_grid(. ~ subtype, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(x = "Drugs", y = "Samples") + 
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(colour = "grey", fill = "white"),
    strip.text.x = element_text(angle = 90),
    strip.text.y.left = element_text(face = 'bold', angle = 0),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = rel(RELSIZE)),
    legend.text = element_text(size = rel(RELSIZE)),
    legend.key.height = unit(1, "cm"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    text = element_text(size = FONTSIZE)
  )

p
my_save_plot(p, 'plots_manuscript/mp_specific_drugs.pdf', width = 9, height = 10)
