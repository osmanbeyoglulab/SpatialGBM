setwd("~/Desktop/Pitt/Project_zhang_gbm_v2/")

library(dplyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)
source('aux_correlation.R')

drugs = list.files(path = 'dataframes/correlation_ravi/drug_with_drug/', full.names = FALSE) %>% 
  str_remove(".csv")
dir.create('plots_website/corr_with_drug')
pname = 'plots_website/corr_with_drug/%s.png'
fname = 'dataframes/correlation_ravi/drug_with_%s/%s.csv'

# plot metaprogram proportion
df_mp = read.csv('dataframes/metaprogram_proportion.csv', row.names=1)

# ========== ========== ========== ========== ==========

for (drug in drugs){
  # read files
  df_cor_lig = my_read_csv(sprintf(fname, 'ligand', drug))
  df_cor_rec = my_read_csv(sprintf(fname, 'receptor', drug))
  df_cor_tf = my_read_csv(sprintf(fname, 'tf', drug))
  df_cor_pw = my_read_csv(sprintf(fname, 'pathway', drug), is_pathway = T)
  df_cor_gene = my_read_csv(sprintf(fname, 'gene', drug))
  df_cor_drug = my_read_csv(sprintf(fname, 'drug', drug))
  rownames(df_cor_drug) = str_extract(rownames(df_cor_drug), "^[^|]+")
  
  ligands2 = filter_rows(df_cor_lig, use_abs = F)
  receptors2 = filter_rows(df_cor_rec, use_abs = F)
  pws2 = filter_rows(df_cor_pw, ntop = 8, use_abs = F)
  tfs2 = filter_rows(df_cor_tf, use_abs = F)
  drugs2 = filter_rows(df_cor_drug, threshold = 1, use_abs = F)
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

