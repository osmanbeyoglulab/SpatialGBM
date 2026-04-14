setwd("~/Desktop/Pitt/Project_zhang_gbm_v2/")

library(dplyr)
library(MatrixGenerics)
library(stringr)
library(circlize)
library(ComplexHeatmap)
source('aux_correlation.R')

tfs = list.files(path = 'dataframes/correlation_ravi/tf_with_tf/', full.names = FALSE) %>% 
  str_remove(".csv") 
fname = 'dataframes/correlation_ravi/tf_with_%s/%s.csv'

# plot metaprogram proportion
df_mp = read.csv('dataframes/metaprogram_proportion.csv', row.names=1)

# ========== ========== ========== ========== ==========

tf = 'HIF1A'
df_cor_lig = my_read_csv(sprintf(fname, 'ligand', tf))
df_cor_rec = my_read_csv(sprintf(fname, 'receptor', tf))
df_cor_pw = my_read_csv(sprintf(fname, 'pathway', tf), is_pathway = T)
df_cor_tf = my_read_csv(sprintf(fname, 'tf', tf))
df_cor_gene = my_read_csv(sprintf(fname, 'gene', tf))
df_cor_drug = my_read_csv(sprintf(fname, 'drug', tf))
rownames(df_cor_drug) = str_extract(rownames(df_cor_drug), "^[^|]+")

ligands2 = filter_rows(df_cor_lig, ntop = 15)
receptors2 = filter_rows(df_cor_rec, ntop = 15)
pws2 = filter_rows(df_cor_pw, ntop = 6)
tfs2 = filter_rows(df_cor_tf, ntop = 15)
drugs2 = filter_rows(df_cor_drug, ntop = 15, use_abs = F)
genes2 = filter_rows(df_cor_gene, ntop = 15,
                     gene2remove=rownames(df_cor_lig) %>% union(rownames(df_cor_rec)))
c(ligands2,receptors2,genes2, tf)

n = length(c(ligands2, receptors2, genes2, tfs2, pws2, drugs2))
samples = colnames(df_cor_lig) %>% intersect(colnames(df_cor_rec)) %>%
  intersect(colnames(df_cor_tf)) %>% intersect(colnames(df_cor_gene)) %>%
  intersect(colnames(df_cor_pw)) %>% intersect(colnames(df_cor_drug))
lim = 0.8

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

p

my_save_pdf(p, 'plots_manuscript/corr_with_tf_HIF1A.pdf', width = 19.5, height = 7.5)

# ========== ========== ========== ========== ==========

## For manuscript
tf = 'SOX2'
df_cor_lig = my_read_csv(sprintf(fname, 'ligand', tf))
df_cor_rec = my_read_csv(sprintf(fname, 'receptor', tf))
df_cor_pw = my_read_csv(sprintf(fname, 'pathway', tf), is_pathway = T)
df_cor_tf = my_read_csv(sprintf(fname, 'tf', tf))
df_cor_gene = my_read_csv(sprintf(fname, 'gene', tf))
df_cor_drug = my_read_csv(sprintf(fname, 'drug', tf))
rownames(df_cor_drug) = str_extract(rownames(df_cor_drug), "^[^|]+")

ligands2 = filter_rows(df_cor_lig, ntop = 15)
receptors2 = filter_rows(df_cor_rec, ntop = 15)
pws2 = filter_rows(df_cor_pw, ntop = 6)
tfs2 = filter_rows(df_cor_tf, ntop = 15)
drugs2 = filter_rows(df_cor_drug, ntop = 15, use_abs = F)
genes2 = filter_rows(df_cor_gene, ntop = 15,
                     gene2remove=rownames(df_cor_lig) %>% union(rownames(df_cor_rec)))
c(ligands2,receptors2,genes2, tf)

n = length(c(ligands2, receptors2, genes2, tfs2, pws2, drugs2))
samples = colnames(df_cor_lig) %>% intersect(colnames(df_cor_rec)) %>%
  intersect(colnames(df_cor_tf)) %>% intersect(colnames(df_cor_gene)) %>%
  intersect(colnames(df_cor_pw)) %>% intersect(colnames(df_cor_drug))
lim = 0.8

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

p

my_save_pdf(p, 'plots_manuscript/corr_with_tf_SOX2.pdf', width = 18.75, height = 7.5)

# ========== ========== ========== ========== ==========

## For manuscript
tf = 'STAT3'
df_cor_lig = my_read_csv(sprintf(fname, 'ligand', tf))
df_cor_rec = my_read_csv(sprintf(fname, 'receptor', tf))
df_cor_pw = my_read_csv(sprintf(fname, 'pathway', tf), is_pathway = T)
df_cor_tf = my_read_csv(sprintf(fname, 'tf', tf))
df_cor_gene = my_read_csv(sprintf(fname, 'gene', tf))
df_cor_drug = my_read_csv(sprintf(fname, 'drug', tf))
rownames(df_cor_drug) = str_extract(rownames(df_cor_drug), "^[^|]+")

ligands2 = filter_rows(df_cor_lig, ntop = 15)
receptors2 = filter_rows(df_cor_rec, ntop = 15)
pws2 = filter_rows(df_cor_pw, ntop = 6)
tfs2 = filter_rows(df_cor_tf, ntop = 15)
drugs2 = filter_rows(df_cor_drug, ntop = 15, use_abs = F)
genes2 = filter_rows(df_cor_gene, ntop = 15,
                     gene2remove=rownames(df_cor_lig) %>% union(rownames(df_cor_rec)))
c(ligands2,receptors2,genes2, tf)

n = length(c(ligands2, receptors2, genes2, tfs2, pws2, drugs2))
samples = colnames(df_cor_lig) %>% intersect(colnames(df_cor_rec)) %>%
  intersect(colnames(df_cor_tf)) %>% intersect(colnames(df_cor_gene)) %>%
  intersect(colnames(df_cor_pw)) %>% intersect(colnames(df_cor_drug))
lim = 0.8

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

p

my_save_pdf(p, 'plots_manuscript/corr_with_tf_STAT3.pdf', width = 15, height = 7.5)

