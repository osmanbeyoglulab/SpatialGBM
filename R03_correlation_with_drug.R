setwd("~/Desktop/Pitt/Project_zhang_gbm_v2/")

library(dplyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)
source('aux_correlation.R')

# ========== ========== ========== ========== ==========

drugs = list.files(path = 'dataframes/correlation_ravi/drug_with_drug/', full.names = FALSE) %>% 
  str_remove(".csv")
dir.create('plots_website/correlation_ravi/')
dir.create('plots_website/correlation_ravi/corr_with_drug')
pname = 'plots_website/correlation_ravi/corr_with_drug/%s.png'
fname = 'dataframes/correlation_ravi/drug_with_%s/%s.csv'

# plot metaprogram proportion
df_mp = read.csv('dataframes/metaprogram_proportion.csv', row.names=1)
threshold = 0.5

for (drug in drugs){
  # read files
  df_cor_lig = my_read_csv(sprintf(fname, 'ligand', drug))
  df_cor_rec = my_read_csv(sprintf(fname, 'receptor', drug))
  df_cor_tf = my_read_csv(sprintf(fname, 'tf', drug))
  df_cor_pw = my_read_csv(sprintf(fname, 'pathway', drug), is_pathway = T)
  df_cor_gene = my_read_csv(sprintf(fname, 'gene', drug))
  df_cor_drug = my_read_csv(sprintf(fname, 'drug', drug))
  rownames(df_cor_drug) = str_extract(rownames(df_cor_drug), "^[^|]+")
  
  ligands2 = filter_rows(df_cor_lig, threshold=threshold, use_abs = F)
  receptors2 = filter_rows(df_cor_rec, threshold=threshold, use_abs = F)
  pws2 = filter_rows(df_cor_pw, ntop = 8, threshold=threshold, use_abs = F)
  tfs2 = filter_rows(df_cor_tf, threshold=threshold, use_abs = F)
  drugs2 = filter_rows(df_cor_drug, threshold=threshold, use_abs = F)
  genes2 = filter_rows(df_cor_gene, 
                       gene2remove=rownames(df_cor_lig) %>% union(rownames(df_cor_rec)), use_abs = F)
  
  n = length(c(ligands2, receptors2, genes2, tfs2, pws2, drugs2))
  samples = colnames(df_cor_lig) %>% intersect(colnames(df_cor_rec)) %>% 
    intersect(colnames(df_cor_tf)) %>% intersect(colnames(df_cor_gene)) %>% 
    intersect(colnames(df_cor_pw)) %>% intersect(colnames(df_cor_drug))
  lim = 0.8
  
  if (n>0){
    p = heatmap_proportion(df_mp[, samples])
    drug_rename = gsub("_", "\n", drug)
    if (length(drugs2)>0){
      p = heatmap_drug(df_cor_drug[drugs2, samples], drug_rename, 'Drugs', lim = 1) + p
    }
    if (length(genes2)>0){
      p = heatmap_by_sample(df_cor_gene[genes2, samples], drug_rename, 'Target Genes', lim = lim) + p
    }
    if (length(tfs2)>0){
      p = heatmap_by_sample(df_cor_tf[tfs2, samples], drug_rename, 'TFs', lim = lim) + p
    }
    if (length(pws2)>0){
      p = heatmap_by_sample(df_cor_pw[pws2, samples], drug_rename, 'Pathways', lim = lim) + p
    }
    if (length(receptors2)>0){
      p = heatmap_by_sample(df_cor_rec[receptors2, samples], drug_rename, 'Receptors', lim = lim) + p
    }
    if (length(ligands2)>0){
      p = heatmap_by_sample(df_cor_lig[ligands2, samples], drug_rename, 'Ligands', lim = lim) + p
    }
    
    my_save_plot(p, sprintf(pname, drug), 
                 width = 8+0.2*n, ## Please tune if necessary.
                 height = 3.2+0.2*length(samples)) ## Please tune if necessary.
  }
}

# ========== ========== ========== ========== ==========

drugs = list.files(path = 'dataframes/correlation_ren/drug_with_drug/', full.names = FALSE) %>% 
  str_remove(".csv")
dir.create('plots_website/correlation_ren/')
dir.create('plots_website/correlation_ren/corr_with_drug')
pname = 'plots_website/correlation_ren/corr_with_drug/%s.png'
fname = 'dataframes/correlation_ren/drug_with_%s/%s.csv'

for (drug in drugs){
  # read files
  df_cor_lig = my_read_csv(sprintf(fname, 'ligand', drug))
  df_cor_rec = my_read_csv(sprintf(fname, 'receptor', drug))
  df_cor_tf = my_read_csv(sprintf(fname, 'tf', drug))
  df_cor_pw = my_read_csv(sprintf(fname, 'pathway', drug), is_pathway = T)
  df_cor_gene = my_read_csv(sprintf(fname, 'gene', drug))
  df_cor_drug = my_read_csv(sprintf(fname, 'drug', drug))
  rownames(df_cor_drug) = str_extract(rownames(df_cor_drug), "^[^|]+")
  
  ligands2 = filter_rows(df_cor_lig, threshold=threshold, use_abs = F, ncount = 4)
  receptors2 = filter_rows(df_cor_rec, threshold=threshold, use_abs = F, ncount = 4)
  pws2 = filter_rows(df_cor_pw, ntop = 8, threshold=threshold, use_abs = F, ncount = 4)
  tfs2 = filter_rows(df_cor_tf, threshold=threshold, use_abs = F, ncount = 4)
  drugs2 = filter_rows(df_cor_drug, threshold=threshold, use_abs = F, ncount = 4)
  genes2 = filter_rows(df_cor_gene, ncount = 4, 
                       gene2remove=rownames(df_cor_lig) %>% union(rownames(df_cor_rec)), use_abs = F)
  
  n = length(c(ligands2, receptors2, genes2, tfs2, pws2, drugs2))
  samples = colnames(df_cor_lig) %>% intersect(colnames(df_cor_rec)) %>% 
    intersect(colnames(df_cor_tf)) %>% intersect(colnames(df_cor_gene)) %>% 
    intersect(colnames(df_cor_pw)) %>% intersect(colnames(df_cor_drug))
  lim = 0.8
  
  if (n>0){
    p <- NULL
    drug_rename = gsub("_", "\n", drug) %>% paste0('\n')
    # if ((length(drugs2) > 0) & (length(drugs2) < n)) {
    if (length(drugs2)>0){
      p <- heatmap_drug(df_cor_drug[drugs2, samples], drug_rename, 'Drugs', lim = 1)
    }
    if (length(genes2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_gene[genes2, samples], drug_rename, 'Target Genes', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_gene[genes2, samples], drug_rename, 'Target Genes', lim = lim) + p
      }
    }
    if (length(tfs2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_tf[tfs2, samples], drug_rename, 'TFs', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_tf[tfs2, samples], drug_rename, 'TFs', lim = lim) + p
      }
    }
    if (length(pws2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_pw[pws2, samples], drug_rename, 'Pathways', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_pw[pws2, samples], drug_rename, 'Pathways', lim = lim) + p
      }
    }
    if (length(receptors2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_rec[receptors2, samples], drug_rename, 'Receptors', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_rec[receptors2, samples], drug_rename, 'Receptors', lim = lim) + p
      }
    }
    if (length(ligands2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_lig[ligands2, samples], drug_rename, 'Ligands', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_lig[ligands2, samples], drug_rename, 'Ligands', lim = lim) + p
      }
    }
    
    my_save_plot(p, sprintf(pname, drug), 
                 width = 6+0.2*n, ## Please tune if necessary.
                 height = 6.5) ## Please tune if necessary.
  }
}

# ========== ========== ========== ========== ==========

drugs = list.files(path = 'dataframes/correlation_ravi/drug_with_drug/', full.names = FALSE) %>% 
  str_remove(".csv") %>% 
  intersect(list.files(path = 'dataframes/correlation_ren/drug_with_drug/', full.names = FALSE) %>% 
              str_remove(".csv"))
dir.create('plots_website/correlation_comb/')
dir.create('plots_website/correlation_comb/corr_with_drug')
pname = 'plots_website/correlation_comb/corr_with_drug/%s.png'

fname_ravi = 'dataframes/correlation_ravi/drug_with_%s/%s.csv'
fname_ren = 'dataframes/correlation_ren/drug_with_%s/%s.csv'
ncount = 9

for (drug in drugs){
  # read files
  df_cor_lig = my_combine_csv(c(sprintf(fname_ravi, 'ligand', drug), sprintf(fname_ren, 'ligand', drug)))
  df_cor_rec = my_combine_csv(c(sprintf(fname_ravi, 'receptor', drug), sprintf(fname_ren, 'receptor', drug)))
  df_cor_gene = my_combine_csv(c(sprintf(fname_ravi, 'gene', drug), sprintf(fname_ren, 'gene', drug)))
  df_cor_tf = my_combine_csv(c(sprintf(fname_ravi, 'tf', drug), sprintf(fname_ren, 'tf', drug)))
  df_cor_pw = my_combine_csv(c(sprintf(fname_ravi, 'pathway', drug), sprintf(fname_ren, 'pathway', drug)), 
                             is_pathway = T)
  df_cor_drug = my_combine_csv(c(sprintf(fname_ravi, 'drug', drug), sprintf(fname_ren, 'drug', drug)))
  rownames(df_cor_drug) = str_extract(rownames(df_cor_drug), "^[^|]+")
  
  ligands2 = filter_rows(df_cor_lig, threshold=threshold, use_abs = F, ncount = ncount)
  receptors2 = filter_rows(df_cor_rec, threshold=threshold, use_abs = F, ncount = ncount)
  pws2 = filter_rows(df_cor_pw, ntop = 8, threshold=threshold, use_abs = F, ncount = ncount)
  tfs2 = filter_rows(df_cor_tf, threshold=threshold, use_abs = F, ncount = ncount)
  drugs2 = filter_rows(df_cor_drug, threshold=threshold, use_abs = F, ncount = ncount)
  genes2 = filter_rows(df_cor_gene, ncount = ncount, 
                       gene2remove=rownames(df_cor_lig) %>% union(rownames(df_cor_rec)), use_abs = F)
  
  n = length(c(ligands2, receptors2, genes2, tfs2, pws2, drugs2))
  samples = colnames(df_cor_lig) %>% intersect(colnames(df_cor_rec)) %>% 
    intersect(colnames(df_cor_tf)) %>% intersect(colnames(df_cor_gene)) %>% 
    intersect(colnames(df_cor_pw)) %>% intersect(colnames(df_cor_drug))
  lim = 0.8
  
  if (n>0){
    p <- NULL
    drug_rename = gsub("_", "\n", drug) %>% paste0('\n')
    # if ((length(drugs2) > 0) & (length(drugs2) < n)) {
    if (length(drugs2)>0){
      p <- heatmap_drug(df_cor_drug[drugs2, samples], drug_rename, 'Drugs', lim = 1)
    }
    if (length(genes2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_gene[genes2, samples], drug_rename, 'Target Genes', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_gene[genes2, samples], drug_rename, 'Target Genes', lim = lim) + p
      }
    }
    if (length(tfs2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_tf[tfs2, samples], drug_rename, 'TFs', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_tf[tfs2, samples], drug_rename, 'TFs', lim = lim) + p
      }
    }
    if (length(pws2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_pw[pws2, samples], drug_rename, 'Pathways', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_pw[pws2, samples], drug_rename, 'Pathways', lim = lim) + p
      }
    }
    if (length(receptors2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_rec[receptors2, samples], drug_rename, 'Receptors', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_rec[receptors2, samples], drug_rename, 'Receptors', lim = lim) + p
      }
    }
    if (length(ligands2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_lig[ligands2, samples], drug_rename, 'Ligands', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_lig[ligands2, samples], drug_rename, 'Ligands', lim = lim) + p
      }
    }
    
    my_save_plot(p, sprintf(pname, drug), 
                 width = 6+0.2*n, ## Please tune if necessary.
                 height = 12) ## Please tune if necessary.
  }
}

