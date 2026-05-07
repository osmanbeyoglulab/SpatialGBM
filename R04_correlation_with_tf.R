setwd("~/Desktop/Pitt/Project_zhang_gbm_v2/")

library(dplyr)
library(MatrixGenerics)
library(stringr)
library(circlize)
library(ComplexHeatmap)
source('aux_correlation.R')

# ========== ========== ========== ========== ==========

tfs = list.files(path = 'dataframes/correlation_ravi/tf_with_tf', full.names = FALSE) %>% 
  str_remove(".csv") 
dir.create('plots_website/correlation_ravi/')
dir.create('plots_website/correlation_ravi/corr_with_tf')
pname = 'plots_website/correlation_ravi/corr_with_tf/%s.png'
fname = 'dataframes/correlation_ravi/tf_with_%s/%s.csv'

# plot metaprogram proportion
df_mp = read.csv('dataframes/metaprogram_proportion.csv', row.names=1)

for (tf in tfs){
  # read files
  df_cor_lig = my_read_csv(sprintf(fname, 'ligand', tf))
  df_cor_rec = my_read_csv(sprintf(fname, 'receptor', tf))
  df_cor_pw = my_read_csv(sprintf(fname, 'pathway', tf), is_pathway = T)
  df_cor_tf = my_read_csv(sprintf(fname, 'tf', tf))
  df_cor_gene = my_read_csv(sprintf(fname, 'gene', tf))
  df_cor_drug = my_read_csv(sprintf(fname, 'drug', tf))
  rownames(df_cor_drug) = str_extract(rownames(df_cor_drug), "^[^|]+")
  
  ligands2 = filter_rows(df_cor_lig)
  receptors2 = filter_rows(df_cor_rec)
  pws2 = filter_rows(df_cor_pw, ntop = 8)
  tfs2 = filter_rows(df_cor_tf)
  drugs2 = filter_rows(df_cor_drug, use_abs = F)
  genes2 = filter_rows(df_cor_gene,
                       gene2remove=rownames(df_cor_lig) %>% union(rownames(df_cor_rec)))
  
  n = length(c(ligands2, receptors2, genes2, tfs2, pws2, drugs2))
  samples = colnames(df_cor_lig) %>% intersect(colnames(df_cor_rec)) %>% 
    intersect(colnames(df_cor_tf)) %>% intersect(colnames(df_cor_gene)) %>% 
    intersect(colnames(df_cor_pw)) %>% intersect(colnames(df_cor_drug))
  lim = 0.8
  
  if (n>0){
    p = heatmap_proportion(df_mp[, samples])
    if (length(drugs2)>0){
      p = heatmap_drug(df_cor_drug[drugs2, samples], tf, 'Drugs', lim = 1) + p
    }
    if (length(genes2)>0){
      p = heatmap_by_sample(df_cor_gene[genes2, samples], tf, 'Target Genes', lim = lim) + p
    }
    if (length(tfs2)>0){
      p = heatmap_by_sample(df_cor_tf[tfs2, samples], tf, 'TFs', lim = lim) + p
    }
    if (length(pws2)>0){
      p = heatmap_by_sample(df_cor_pw[pws2, samples], tf, 'Pathways', lim = lim) + p
    }
    if (length(receptors2)>0){
      p = heatmap_by_sample(df_cor_rec[receptors2, samples], tf, 'Receptors', lim = lim) + p
    }
    if (length(ligands2)>0){
      p = heatmap_by_sample(df_cor_lig[ligands2, samples], tf, 'Ligands', lim = lim) + p
    }
    
    my_save_plot(p, sprintf(pname, tf), 
                 width = 8+0.2*n, ## Please tune if necessary.
                 height = 8) ## Please tune if necessary.
  }
}

# ========== ========== ========== ========== ==========

tfs = list.files(path = 'dataframes/correlation_ren/tf_with_tf', full.names = FALSE) %>% 
  str_remove(".csv") 
dir.create('plots_website/correlation_ren/')
dir.create('plots_website/correlation_ren/corr_with_tf')
pname = 'plots_website/correlation_ren/corr_with_tf/%s.png'
fname = 'dataframes/correlation_ren/tf_with_%s/%s.csv'

for (tf in tfs){
  # read files
  df_cor_lig = my_read_csv(sprintf(fname, 'ligand', tf))
  df_cor_rec = my_read_csv(sprintf(fname, 'receptor', tf))
  df_cor_pw = my_read_csv(sprintf(fname, 'pathway', tf), is_pathway = T)
  df_cor_tf = my_read_csv(sprintf(fname, 'tf', tf))
  df_cor_gene = my_read_csv(sprintf(fname, 'gene', tf))
  df_cor_drug = my_read_csv(sprintf(fname, 'drug', tf))
  rownames(df_cor_drug) = str_extract(rownames(df_cor_drug), "^[^|]+")
  
  ligands2 = filter_rows(df_cor_lig, ncount = 4)
  receptors2 = filter_rows(df_cor_rec, ncount = 4)
  pws2 = filter_rows(df_cor_pw, ntop = 8, ncount = 4)
  tfs2 = filter_rows(df_cor_tf, ncount = 4)
  drugs2 = filter_rows(df_cor_drug, use_abs = F, ncount = 4)
  genes2 = filter_rows(df_cor_gene, ncount = 4,
                       gene2remove=rownames(df_cor_lig) %>% union(rownames(df_cor_rec)))
  
  n = length(c(ligands2, receptors2, genes2, tfs2, pws2, drugs2))
  samples = colnames(df_cor_lig) %>% intersect(colnames(df_cor_rec)) %>% 
    intersect(colnames(df_cor_tf)) %>% intersect(colnames(df_cor_gene)) %>% 
    intersect(colnames(df_cor_pw)) %>% intersect(colnames(df_cor_drug))
  lim = 0.8
  
  if (n>0){
    p <- NULL
    if (length(drugs2) > 0) {
      p <- heatmap_drug(df_cor_drug[drugs2, samples], tf, 'Drugs', lim = 1)
    }
    if (length(genes2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_gene[genes2, samples], tf, 'Target Genes', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_gene[genes2, samples], tf, 'Target Genes', lim = lim) + p
      }
    }
    if (length(tfs2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_tf[tfs2, samples], tf, 'TFs', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_tf[tfs2, samples], tf, 'TFs', lim = lim) + p
      }
    }
    if (length(pws2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_pw[pws2, samples], tf, 'Pathways', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_pw[pws2, samples], tf, 'Pathways', lim = lim) + p
      }
    }
    if (length(receptors2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_rec[receptors2, samples], tf, 'Receptors', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_rec[receptors2, samples], tf, 'Receptors', lim = lim) + p
      }
    }
    if (length(ligands2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_lig[ligands2, samples], tf, 'Ligands', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_lig[ligands2, samples], tf, 'Ligands', lim = lim) + p
      }
    }
    
    my_save_plot(p, sprintf(pname, tf), 
                 width = 5+0.2*n, ## Please tune if necessary.
                 height = 6.5) ## Please tune if necessary.
  }
}

# ========== ========== ========== ========== ==========

tfs = list.files(path = 'dataframes/correlation_ravi/tf_with_tf/', full.names = FALSE) %>% 
  str_remove(".csv") %>% 
  intersect(list.files(path = 'dataframes/correlation_ren/tf_with_tf/', full.names = FALSE) %>% 
              str_remove(".csv"))
dir.create('plots_website/correlation_comb/')
dir.create('plots_website/correlation_comb/corr_with_tf')
pname = 'plots_website/correlation_comb/corr_with_tf/%s.png'

fname_ravi = 'dataframes/correlation_ravi/tf_with_%s/%s.csv'
fname_ren = 'dataframes/correlation_ren/tf_with_%s/%s.csv'
ncount = 9

for (tf in tfs){
  # read files
  df_cor_lig = my_combine_csv(c(sprintf(fname_ravi, 'ligand', tf), sprintf(fname_ren, 'ligand', tf)))
  df_cor_rec = my_combine_csv(c(sprintf(fname_ravi, 'receptor', tf), sprintf(fname_ren, 'receptor', tf)))
  df_cor_gene = my_combine_csv(c(sprintf(fname_ravi, 'gene', tf), sprintf(fname_ren, 'gene', tf)))
  df_cor_tf = my_combine_csv(c(sprintf(fname_ravi, 'tf', tf), sprintf(fname_ren, 'tf', tf)))
  df_cor_pw = my_combine_csv(c(sprintf(fname_ravi, 'pathway', tf), sprintf(fname_ren, 'pathway', tf)), 
                             is_pathway = T)
  df_cor_drug = my_combine_csv(c(sprintf(fname_ravi, 'drug', tf), sprintf(fname_ren, 'drug', tf)))
  rownames(df_cor_drug) = str_extract(rownames(df_cor_drug), "^[^|]+")
  
  ligands2 = filter_rows(df_cor_lig, ncount = ncount)
  receptors2 = filter_rows(df_cor_rec, ncount = ncount)
  pws2 = filter_rows(df_cor_pw, ntop = 8, ncount = ncount)
  tfs2 = filter_rows(df_cor_tf, ncount = ncount)
  drugs2 = filter_rows(df_cor_drug, use_abs = F, ncount = ncount)
  genes2 = filter_rows(df_cor_gene, ncount = ncount,
                       gene2remove=rownames(df_cor_lig) %>% union(rownames(df_cor_rec)))
  
  n = length(c(ligands2, receptors2, genes2, tfs2, pws2, drugs2))
  samples = colnames(df_cor_lig) %>% intersect(colnames(df_cor_rec)) %>% 
    intersect(colnames(df_cor_tf)) %>% intersect(colnames(df_cor_gene)) %>% 
    intersect(colnames(df_cor_pw)) %>% intersect(colnames(df_cor_drug))
  lim = 0.8
  
  if (n>0){
    p <- NULL
    if (length(drugs2) > 0) {
      p <- heatmap_drug(df_cor_drug[drugs2, samples], tf, 'Drugs', lim = 1)
    }
    if (length(genes2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_gene[genes2, samples], tf, 'Target Genes', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_gene[genes2, samples], tf, 'Target Genes', lim = lim) + p
      }
    }
    if (length(tfs2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_tf[tfs2, samples], tf, 'TFs', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_tf[tfs2, samples], tf, 'TFs', lim = lim) + p
      }
    }
    if (length(pws2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_pw[pws2, samples], tf, 'Pathways', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_pw[pws2, samples], tf, 'Pathways', lim = lim) + p
      }
    }
    if (length(receptors2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_rec[receptors2, samples], tf, 'Receptors', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_rec[receptors2, samples], tf, 'Receptors', lim = lim) + p
      }
    }
    if (length(ligands2) > 0) {
      if (is.null(p)) {
        p <- heatmap_by_sample(df_cor_lig[ligands2, samples], tf, 'Ligands', lim = lim)
      } else {
        p <- heatmap_by_sample(df_cor_lig[ligands2, samples], tf, 'Ligands', lim = lim) + p
      }
    }
    
    my_save_plot(p, sprintf(pname, tf), 
                 width = 5+0.2*n, ## Please tune if necessary.
                 height = 12) ## Please tune if necessary.
  }
}
