setwd("~/Desktop/Pitt/Project_zhang_gbm_v2/")

library(dplyr)
library(MatrixGenerics)
library(stringr)
library(circlize)
library(ComplexHeatmap)
source('aux_correlation.R')

# ========== ========== ========== ========== ==========

pathways = list.files(path = 'dataframes/correlation_ravi/pathway_with_pathway/', full.names = FALSE) %>% 
  str_remove(".csv")
dir.create('plots_website/correlation_ravi/')
dir.create('plots_website/correlation_ravi/corr_with_pathway')
pname = 'plots_website/correlation_ravi/corr_with_pathway/%s.png'
fname = 'dataframes/correlation_ravi/pathway_with_%s/%s.csv'

# plot metaprogram proportion
df_mp = read.csv('dataframes/metaprogram_proportion.csv', row.names=1)

for (pathway in pathways){
  # read files
  df_cor_lig = my_read_csv(sprintf(fname, 'ligand', pathway))
  df_cor_rec = my_read_csv(sprintf(fname, 'receptor', pathway))
  df_cor_tf = my_read_csv(sprintf(fname, 'tf', pathway))
  df_cor_pw = my_read_csv(sprintf(fname, 'pathway', pathway), is_pathway = T)
  df_cor_gene = my_read_csv(sprintf(fname, 'gene', pathway))
  df_cor_drug = my_read_csv(sprintf(fname, 'drug', pathway))
  rownames(df_cor_drug) = str_extract(rownames(df_cor_drug), "^[^|]+")
  
  ligands2 = filter_rows(df_cor_lig)
  receptors2 = filter_rows(df_cor_rec)
  tfs2 = filter_rows(df_cor_tf)
  pws2 = filter_rows(df_cor_pw, ntop = 8)
  drugs2 = filter_rows(df_cor_drug, ntop = 8, use_abs = F)
  genes2 = filter_rows(df_cor_gene,
                       gene2remove=rownames(df_cor_lig) %>% union(rownames(df_cor_rec)))
  
  n = length(c(ligands2, receptors2, genes2, tfs2, pws2, drugs2))
  samples = colnames(df_cor_lig) %>% intersect(colnames(df_cor_rec)) %>% 
    intersect(colnames(df_cor_tf)) %>% intersect(colnames(df_cor_gene)) %>% 
    intersect(colnames(df_cor_pw)) %>% intersect(colnames(df_cor_drug))
  lim = 1
  pathway2 = hallmarks[pathway]
  
  if (n>0){
    p = heatmap_proportion(df_mp[, samples])
    if (length(drugs2)>0){
      p = heatmap_drug(df_cor_drug[drugs2, samples], pathway2, 'Drugs', lim = lim) + p
    }
    if (length(genes2)>0){
      p = heatmap_by_sample(df_cor_gene[genes2, samples], pathway2, 'Target Genes', lim = lim) + p
    }
    if (length(tfs2)>0){
      p = heatmap_by_sample(df_cor_tf[tfs2, samples], pathway2, 'TFs', lim = lim) + p
    }
    if (length(pws2)>0){
      p = heatmap_by_sample(df_cor_pw[pws2, samples], pathway2, 'Pathways', lim = lim) + p
    }
    if (length(receptors2)>0){
      p = heatmap_by_sample(df_cor_rec[receptors2, samples], pathway2, 'Receptors', lim = lim) + p
    }
    if (length(ligands2)>0){
      p = heatmap_by_sample(df_cor_lig[ligands2, samples], pathway2, 'Ligands', lim = lim) + p
    }
    
    my_save_plot(p, sprintf(pname, pathway), 
                 width = 8+0.2*n, ## Please tune if necessary.
                 height = 3.2+0.2*length(samples)) ## Please tune if necessary.
  }
}

# ========== ========== ========== ========== ==========

pathways = list.files(path = 'dataframes/correlation_ren/pathway_with_pathway/', full.names = FALSE) %>% 
  str_remove(".csv")
dir.create('plots_website/correlation_ren/')
dir.create('plots_website/correlation_ren/corr_with_pathway')
pname = 'plots_website/correlation_ren/corr_with_pathway/%s.png'
fname = 'dataframes/correlation_ren/pathway_with_%s/%s.csv'

for (pathway in pathways){
  # read files
  df_cor_lig = my_read_csv(sprintf(fname, 'ligand', pathway))
  df_cor_rec = my_read_csv(sprintf(fname, 'receptor', pathway))
  df_cor_tf = my_read_csv(sprintf(fname, 'tf', pathway))
  df_cor_pw = my_read_csv(sprintf(fname, 'pathway', pathway), is_pathway = T)
  df_cor_gene = my_read_csv(sprintf(fname, 'gene', pathway))
  df_cor_drug = my_read_csv(sprintf(fname, 'drug', pathway))
  rownames(df_cor_drug) = str_extract(rownames(df_cor_drug), "^[^|]+")
  
  ligands2 = filter_rows(df_cor_lig, ncount = 4)
  receptors2 = filter_rows(df_cor_rec, ncount = 4)
  tfs2 = filter_rows(df_cor_tf, ncount = 4)
  pws2 = filter_rows(df_cor_pw, ntop = 8, ncount = 4)
  drugs2 = filter_rows(df_cor_drug, ntop = 8, use_abs = F, ncount = 4)
  genes2 = filter_rows(df_cor_gene, ncount = 4,
                       gene2remove=rownames(df_cor_lig) %>% union(rownames(df_cor_rec)))
  
  n = length(c(ligands2, receptors2, genes2, tfs2, pws2, drugs2))
  samples = colnames(df_cor_lig) %>% intersect(colnames(df_cor_rec)) %>% 
    intersect(colnames(df_cor_tf)) %>% intersect(colnames(df_cor_gene)) %>% 
    intersect(colnames(df_cor_pw)) %>% intersect(colnames(df_cor_drug))
  lim = 1
  pathway2 = hallmarks[pathway]
  
  if (n>0){
    p <- NULL
    if (length(drugs2) > 0) {
      p <- heatmap_drug(df_cor_drug[drugs2, samples], pathway2, 'Drugs', lim = 1)
    }
    if (length(genes2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_gene[genes2, samples], pathway2, 'Target Genes', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_gene[genes2, samples], pathway2, 'Target Genes', lim = lim) + p
      }
    }
    if (length(tfs2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_tf[tfs2, samples], pathway2, 'TFs', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_tf[tfs2, samples], pathway2, 'TFs', lim = lim) + p
      }
    }
    if (length(pws2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_pw[pws2, samples], pathway2, 'Pathways', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_pw[pws2, samples], pathway2, 'Pathways', lim = lim) + p
      }
    }
    if (length(receptors2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_rec[receptors2, samples], pathway2, 'Receptors', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_rec[receptors2, samples], pathway2, 'Receptors', lim = lim) + p
      }
    }
    if (length(ligands2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_lig[ligands2, samples], pathway2, 'Ligands', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_lig[ligands2, samples], pathway2, 'Ligands', lim = lim) + p
      }
    }
    
    my_save_plot(p, sprintf(pname, pathway), 
                 width = 5+0.2*n, ## Please tune if necessary.
                 height = 6.5) ## Please tune if necessary.
  }
}
