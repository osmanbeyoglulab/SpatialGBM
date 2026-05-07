cmap = list(
  'Ligands' = c("#2166AC", "white", "#B2182B"), # brewer.pal(n=11, name = "RdBu")
  'Receptors' = c("#4D9221", "white", "#C51B7D"), # brewer.pal(n=11, name = "PiYG")
  'Pathways' = c("#542788", "white", "#B35806"), # brewer.pal(n=11, name = "PuOr")
  'TFs' = c("#01665E", "white", "#8C510A"), # brewer.pal(n=11, name = "BrBG")
  'Target Genes' = c("#4575B4", "white", "#D73027") # brewer.pal(n=11, name = "RdYlBu")
  # 'Drugs' = c("#4292C6", "white", "#F16913") # "Blues" and "Oranges"
)
# [1] "#FFFFCC" "#FFEDA0" "#FED976" "#FEB24C" "#FD8D3C" "#FC4E2A" "#E31A1C" "#BD0026"
# [9] "#800026"

metaprogram_palette <- c(
  'AC' = '#FF7F00',
  'Chromatin.Reg' = '#3969AC', 
  'Inflammatory.Mac' = '#11A579', 
  'MES' = '#F6CF71', 
  'MES.Ast' = '#66C5CC', 
  'MES.Hyp' = '#F2B701',
  'Mac' = '#80BA5A', 
  'NPC' = '#7F3C8D', 
  'Neuron' = '#4b4b8f', 
  'OPC' = '#f97b72', 
  'Oligo' = '#ACA4E2', 
  'Prolif.Metab' = '#E73F74', 
  'Reactive.Ast' = '#D4D915',
  'Vasc' = '#CF1C90'
)

hallmarks = list(
  'HALLMARK_ADIPOGENESIS' = 'Adipogenesis', 
  'HALLMARK_ALLOGRAFT_REJECTION' = 'Allograft Rejection',
  'HALLMARK_ANDROGEN_RESPONSE' = 'Androgen Response', 
  'HALLMARK_ANGIOGENESIS' = 'Angiogenesis', 
  'HALLMARK_APICAL_JUNCTION' = 'Apical Junction', 
  'HALLMARK_APICAL_SURFACE' = 'Apical Surface', 
  'HALLMARK_APOPTOSIS' = 'Apoptosis', 
  'HALLMARK_BILE_ACID_METABOLISM' = 'Bile Acid Metabolism',
  'HALLMARK_CHOLESTEROL_HOMEOSTASIS' = 'Cholesterol Homeostasis', 
  'HALLMARK_COAGULATION' = 'Coagulation',
  'HALLMARK_COMPLEMENT' = 'Complement',
  'HALLMARK_DNA_REPAIR' = 'DNA Repair', 
  'HALLMARK_E2F_TARGETS' = 'E2F Targets', 
  'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION' = 'EMT', 
  'HALLMARK_ESTROGEN_RESPONSE_EARLY' = 'ER Early', 
  'HALLMARK_ESTROGEN_RESPONSE_LATE' = 'ER Late', 
  'HALLMARK_FATTY_ACID_METABOLISM' = 'Fatty Acid Metabolism', 
  'HALLMARK_G2M_CHECKPOINT' = 'G2M Checkpoint',
  'HALLMARK_GLYCOLYSIS' = 'Glycolysis', 
  'HALLMARK_HEDGEHOG_SIGNALING' = 'Hedgehog Signalling',
  'HALLMARK_HEME_METABOLISM' = 'Heme Metabolism', 
  'HALLMARK_HYPOXIA' = 'Hypoxia',
  'HALLMARK_IL2_STAT5_SIGNALING' = 'IL2-STAT5 Signaling',
  'HALLMARK_IL6_JAK_STAT3_SIGNALING' = 'IL6/JAK/STAT3 Signaling',
  'HALLMARK_INFLAMMATORY_RESPONSE' = 'Inflammatory Response', 
  'HALLMARK_INTERFERON_ALPHA_RESPONSE' = 'INFA Response', 
  'HALLMARK_INTERFERON_GAMMA_RESPONSE' = 'INFG Response',
  'HALLMARK_KRAS_SIGNALING_DN' = 'KRAS Signaling Down',
  'HALLMARK_KRAS_SIGNALING_UP' = 'KRAS Signaling Up', 
  'HALLMARK_MITOTIC_SPINDLE' = 'Mitotic Spindle', 
  'HALLMARK_MTORC1_SIGNALING' = 'mTORC1 Signaling',
  'HALLMARK_MYC_TARGETS_V1' = 'MYC Targets v1',
  'HALLMARK_MYC_TARGETS_V2' = 'MYC Targets v2',
  'HALLMARK_MYOGENESIS' = 'Myogenesis', 
  'HALLMARK_NOTCH_SIGNALING' = 'Notch Signaling', 
  'HALLMARK_OXIDATIVE_PHOSPHORYLATION' = 'Oxidative phosphorylation', 
  'HALLMARK_P53_PATHWAY' = 'p53', 
  'HALLMARK_PANCREAS_BETA_CELLS' = 'Pancreas Beta Cells', 
  'HALLMARK_PEROXISOME' = 'Peroxisome', 
  'HALLMARK_PI3K_AKT_MTOR_SIGNALING' = 'PI3K/AKT/mTOR Signaling', 
  'HALLMARK_PROTEIN_SECRETION' = 'Protein Secretion', 
  'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY' = 'ROS',  
  'HALLMARK_SPERMATOGENESIS' = 'Spermatogenesis',
  'HALLMARK_TGF_BETA_SIGNALING' = 'TGFB Signaling',
  'HALLMARK_TNFA_SIGNALING_VIA_NFKB' = 'TNFA Signaling via NFKB',
  'HALLMARK_UNFOLDED_PROTEIN_RESPONSE' = 'UPR',
  'HALLMARK_UV_RESPONSE_DN' = 'UV Response Down',
  'HALLMARK_UV_RESPONSE_UP' = 'UV Response Up', 
  'HALLMARK_WNT_BETA_CATENIN_SIGNALING' = 'Wnt/beta-catenin Signaling', 
  'HALLMARK_XENOBIOTIC_METABOLISM' = 'Xenobiotic Metabolism')

marker_genes <- list(
  "Neuron" = c("SNAP25", "SYT1", "SLC17A7"),
  "Vasc" = c("VWF", "PDGFRB", "ESAM",
             'PECAM1', 'KDR', 
             'RGS5', 'MCAM',
             'COL1A1', 'COL3A1', 'DCN', 'LUM'),
  "Mac" = c("CD14", "CD74", "FCGR3A"),
  "Oligo" = c("PLP1", "MAG", "CLDND1"),
  "Inflammatory-Mac" = c("CCL2", "IL1B", "CCL3", "STAT1"),
  
  "MES" = c("EFEMP1", "CHI3L1", "CD44",
            'ANXA1', 'LGALS3', 'FN1', 'VIM'),
  "MES-Ast" = c("GAP43", "KCNF1", "HOXA9", "FABP7"),
  "MES-Hyp" = c("VEGFA", "NDRG1", "PGK1",
                'CA9', 'CA12', 'LDHA', 'ENO1', 'ADM'),
  "AC" = c("S100B", "AQP1", "TTYH1",
           'GFAP', 'ALDH1L1', 'SLC1A3', 'CLU'),
  "OPC" = c("SOX10", "BCAN", "BCAS1",
            'PDGFRA', 'OLIG1', 'OLIG2', 'CSPG4'),
  "Chromatin-Reg" = c("ATFX", "NFX1", "CCNI", "KMT2E", "SOX4"),
  "Prolif-Metab" = c("CCND1", "SLC16A1", "CTNNB1", "SOX9",
                     'MKI67', 'TOP2A', 'CDK1', 'HMGB2'),
  "Reactive-Ast" = c("AQP4", "GJA1", "GSN", "HLA-DRA"),
  "NPC" = c("INA", "NELL2", "DCX", "STMN2",
            'SOX2', 'NES', 'ASCL1', 'DLL3'),
  
  # "Low-quality" = c("MALAT1"),
  "Microglia" = c('P2RY12', 'CX3CR1', 'TMEM119', 'AIF1'),
  "Macrophages" = c('CD68', 'C1QC', 'CTSD', 'LYZ'),
  "T" =	c('CD3D', 'CD3E',	
          'CD8A', 'CD8B', 'GZMB', 'PRF1',
          'FOXP3', 'IL2RA'),
  "NK" = c('NKG7', 'KLRD1', 'GNLY')
)


sample_list = c('MGH258', 'UKF243', 'UKF248', 'UKF251', 'UKF255', 'UKF259', 'UKF260',
                'UKF266', 'UKF269', 'UKF275', 'UKF296', 'UKF304', 'UKF313', 'UKF334',
                'ZH1007inf', 'ZH1007nec', 'ZH1019T1', 'ZH1019inf', 'ZH8811Abulk', 'ZH8811Bbulk',
                'ZH8812bulk', 'ZH881T1', 'ZH881inf', 'ZH916T1', 'ZH916bulk', 'ZH916inf') %>% sort()

my_read_csv = function(fname, is_pathway=F){
  df = read.csv(fname, row.names=1)
  df[is.na(df)] = 0
  if (is_pathway){
    rownames(df) = hallmarks[rownames(df)]
  }
  return(df)
}

my_combine_csv = function(fname_list, is_pathway=F){
  df_list = list()
  for (i in seq_along(fname_list)){
    df_list[[i]] = my_read_csv(fname_list[i], is_pathway)
  }
  common_rows = Reduce(intersect, lapply(df_list, rownames))
  df_combined = do.call(cbind, lapply(df_list, function(df){
    df[common_rows, , drop=FALSE]
  }))
  return(df_combined)
}

filter_rows = function(df, ntop=20, threshold=0.5, ncount=5, use_abs=T, gene2remove=''){ ### Please tune ###
  # ntop=20: Keep top 20 rows by default (adjustable).
  # threshold=0.5: Minimum absolute value to consider a column "significant."
  # ncount=3: Minimum number of columns exceeding threshold.
  if (dim(df)[2] > 1) {
    ntop = min(ntop, dim(df)[1])
    row_ind = rownames(df) %>% setdiff(gene2remove)
    row_means = rowMeans(df[row_ind, ])
    if (use_abs) {row_means = abs(row_means)}
    # Display a maximum of ntop rows given the threshold
    row_ind1 = row_means %>% sort(decreasing = T) %>% .[1:ntop] %>% names()
    row_ind2 = row_ind1[rowSums(abs(df[row_ind1,])>=threshold)>=ncount]
    return(row_ind2)}
  else {
    return( c() )
  }
}


my_save_pdf = function(p, fname, width=10, height=6){
  pdf(fname, width=width, height=height)
  print(p)
  dev.off()
}

my_save_plot = function(p, fname, width=10, height=6){
  png(fname, width=width, height=height, units="in", res=200)
  print(p)
  dev.off()
}

heatmap_proportion = function(df){
  p = Heatmap(df %>% t(),
              col = colorRamp2(seq(0, 0.45, length.out = 9),
                               RColorBrewer::brewer.pal(n=9, name = "Oranges")),
              cluster_columns = F,
              cluster_rows = T,
              column_names_side ="top",
              row_names_side ="left",
              row_title = 'Sapmles',
              column_title = 'Metaprograms',
              name = "Proportion")
  return(p)
}

heatmap_by_sample = function(df, target, cols, lim=1){
  col_names = sprintf('%s-%s\nSpearman Rho', target, str_sub(cols, end = -2))
  p = Heatmap(df %>% t(),
              col = colorRamp2(c(-lim, 0, lim), cmap[[cols]]),
              cluster_columns = T,
              cluster_rows = T,
              column_names_side ="top",
              row_names_side ="left",
              row_title = 'Samples',
              # column_title = str_sub(cols, end = 2),
              column_title = cols,
              name = col_names)
  return(p)
}

heatmap_drug = function(df, target, cols, lim=1){
  col_names = sprintf('%s-%s\nSpearman Rho', target, str_sub(cols, end = -2))
  p = Heatmap(df %>% t(),
              col = colorRamp2(seq(0, lim, length.out = 9),
                               RColorBrewer::brewer.pal(n=9, name = "YlOrRd")),
              cluster_columns = T,
              cluster_rows = T,
              column_names_side ="top",
              row_names_side ="left",
              row_title = 'Samples',
              # column_title = str_sub(cols, end = 2),
              column_title = cols,
              name = col_names)
  return(p)
}


# filter_dataframe_by_threshold <- function(
#     df,
#     thre_val = 2,
#     thre_sam = 3,
#     na.rm = TRUE
# ) {
#   pass_count <- apply(df, 1, function(x) {
#     sum(x > thre_val, na.rm = na.rm)
#   })
#   filtered_df <- df[pass_count >= thre_sam, , drop = FALSE]
#   return(filtered_df)
# }

filter_dataframe_by_threshold <- function(
    df,
    thre_val = 2,
    thre_sam = 3,
    thre_na = 3,  
    na.rm = TRUE
) {
  if (!is.null(thre_na)) {
    na_count <- apply(df, 1, function(x) {
      sum(is.na(x))
    })
    na_pass <- na_count <= thre_na
  } else {
    na_pass <- rep(TRUE, nrow(df))
  }
  
  pass_count <- apply(df, 1, function(x) {
    sum(x > thre_val, na.rm = na.rm)
  })
  pass_both <- na_pass & (pass_count >= thre_sam)
  filtered_df <- df[pass_both, , drop = FALSE]
  
  return(filtered_df)
}
