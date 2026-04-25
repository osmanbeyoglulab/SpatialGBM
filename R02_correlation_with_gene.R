setwd("~/Desktop/Pitt/Project_zhang_gbm_v2/")

library(dplyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)
source('aux_correlation.R')

# ========== ========== ========== ========== ==========

genes = list.files(path = 'dataframes/correlation_ravi/gene_with_gene/', full.names = FALSE) %>% 
  str_remove(".csv")
dir.create('plots_website/correlation_ravi/')
dir.create('plots_website/correlation_ravi/corr_with_gene')
pname = 'plots_website/correlation_ravi/corr_with_gene/%s.png'
fname = 'dataframes/correlation_ravi/gene_with_%s/%s.csv'

for (gene in genes){
  # read files
  df_cor_lig = my_read_csv(sprintf(fname, 'ligand', gene))
  df_cor_rec = my_read_csv(sprintf(fname, 'receptor', gene))
  df_cor_gene = my_read_csv(sprintf(fname, 'gene', gene))
  df_cor_tf = my_read_csv(sprintf(fname, 'tf', gene))
  df_cor_pw = my_read_csv(sprintf(fname, 'pathway', gene), is_pathway = T)
  df_cor_drug = my_read_csv(sprintf(fname, 'drug', gene))
  rownames(df_cor_drug) = str_extract(rownames(df_cor_drug), "^[^|]+")
  
  ligands2 = filter_rows(df_cor_lig, gene2remove = gene)
  receptors2 = filter_rows(df_cor_rec, gene2remove = gene)
  genes2 = filter_rows(df_cor_gene)
  pws2 = filter_rows(df_cor_pw, ntop = 8)
  tfs2 = filter_rows(df_cor_tf)
  drugs2 = filter_rows(df_cor_drug, use_abs = F)
  
  n = length(c(ligands2, receptors2, genes2, tfs2, pws2, drugs2))
  samples = colnames(df_cor_lig) %>% intersect(colnames(df_cor_rec)) %>% 
    intersect(colnames(df_cor_tf)) %>% intersect(colnames(df_cor_gene)) %>% 
    intersect(colnames(df_cor_pw)) %>% intersect(colnames(df_cor_drug))
  lim = 0.8
  
  if (n>0){
    p = heatmap_proportion(df_mp[, samples])
    if (length(drugs2)>0){
      p = heatmap_drug(df_cor_drug[drugs2, samples], gene, 'Drugs', lim = lim) + p
    }
    if (length(genes2)>0){
      p = heatmap_by_sample(df_cor_gene[genes2, samples], gene, 'Target Genes', lim = lim) + p
    }
    if (length(tfs2)>0){
      p = heatmap_by_sample(df_cor_tf[tfs2, samples], gene, 'TFs', lim = lim) + p
    }
    if (length(pws2)>0){
      p = heatmap_by_sample(df_cor_pw[pws2, samples], gene, 'Pathways', lim = lim) + p
    }
    if (length(receptors2)>0){
      p = heatmap_by_sample(df_cor_rec[receptors2, samples], gene, 'Receptors', lim = lim) + p
    }
    if (length(ligands2)>0){
      p = heatmap_by_sample(df_cor_lig[ligands2, samples], gene, 'Ligands', lim = lim) + p
    }
    
    my_save_plot(p, sprintf(pname, gene), 
                 width = 8+0.2*n, ## Please tune if necessary.
                 height = 3.2+0.2*length(samples)) ## Please tune if necessary.
  }
}

# ========== ========== ========== ========== ==========

genes = list.files(path = 'dataframes/correlation_ren/gene_with_gene/', full.names = FALSE) %>% 
  str_remove(".csv")
dir.create('plots_website/correlation_ren/')
dir.create('plots_website/correlation_ren/corr_with_gene')
pname = 'plots_website/correlation_ren/corr_with_gene/%s.png'
fname = 'dataframes/correlation_ren/gene_with_%s/%s.csv'

for (gene in genes){
  # read files
  df_cor_lig = my_read_csv(sprintf(fname, 'ligand', gene))
  df_cor_rec = my_read_csv(sprintf(fname, 'receptor', gene))
  df_cor_gene = my_read_csv(sprintf(fname, 'gene', gene))
  df_cor_tf = my_read_csv(sprintf(fname, 'tf', gene))
  df_cor_pw = my_read_csv(sprintf(fname, 'pathway', gene), is_pathway = T)
  df_cor_drug = my_read_csv(sprintf(fname, 'drug', gene))
  rownames(df_cor_drug) = str_extract(rownames(df_cor_drug), "^[^|]+")
  
  ligands2 = filter_rows(df_cor_lig, gene2remove = gene, ncount = 4)
  receptors2 = filter_rows(df_cor_rec, gene2remove = gene, ncount = 4)
  genes2 = filter_rows(df_cor_gene, ncount = 4)
  pws2 = filter_rows(df_cor_pw, ntop = 8, ncount = 4)
  tfs2 = filter_rows(df_cor_tf, ncount = 4)
  drugs2 = filter_rows(df_cor_drug, use_abs = F, ncount = 4)
  
  n = length(c(ligands2, receptors2, genes2, tfs2, pws2, drugs2))
  samples = colnames(df_cor_lig) %>% intersect(colnames(df_cor_rec)) %>% 
    intersect(colnames(df_cor_tf)) %>% intersect(colnames(df_cor_gene)) %>% 
    intersect(colnames(df_cor_pw)) %>% intersect(colnames(df_cor_drug))
  lim = 0.8
  
  if (n>0){
    p <- NULL
    if (length(drugs2) > 0) {
      p <- heatmap_drug(df_cor_drug[drugs2, samples], gene, 'Drugs', lim = 1)
    }
    if (length(genes2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_gene[genes2, samples], gene, 'Target Genes', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_gene[genes2, samples], gene, 'Target Genes', lim = lim) + p
      }
    }
    if (length(tfs2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_tf[tfs2, samples], gene, 'TFs', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_tf[tfs2, samples], gene, 'TFs', lim = lim) + p
      }
    }
    if (length(pws2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_pw[pws2, samples], gene, 'Pathways', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_pw[pws2, samples], gene, 'Pathways', lim = lim) + p
      }
    }
    if (length(receptors2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_rec[receptors2, samples], gene, 'Receptors', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_rec[receptors2, samples], gene, 'Receptors', lim = lim) + p
      }
    }
    if (length(ligands2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_lig[ligands2, samples], gene, 'Ligands', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_lig[ligands2, samples], gene, 'Ligands', lim = lim) + p
      }
    }
    
    my_save_plot(p, sprintf(pname, gene), 
                 width = 5+0.2*n, ## Please tune if necessary.
                 height = 6.5) ## Please tune if necessary.
  }
}
