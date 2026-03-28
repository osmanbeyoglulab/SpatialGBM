import numpy as np

# figure related variables
figsize = 2
fontsize = 6
dpi = 300 ## change dpi

# metaprogram related variables
metaprograms = {'AC': 'Other', 
                'Chromatin.Reg': 'Other', 
                'Inflammatory.Mac': 'Other', 
                'Mac': 'Other', 
                'MES': 'Other', 
                'MES.Ast': 'Other', 
                'MES.Hyp': 'Other',
                'NPC': 'Other', 
                'Neuron': 'Other', 
                'OPC': 'Other', 
                'Oligo': 'Other', 
                'Prolif.Metab': 'Other', 
                'Reactive.Ast': 'Other',
                'Vasc': 'Other'}

metaprogram_palette = {'AC': '#FF7F00',
                       'Chromatin.Reg': '#3969AC', 
                       'Inflammatory.Mac': '#11A579', 
                       'MES': '#F6CF71', 
                       'MES.Ast': '#66C5CC', 
                       'MES.Hyp': '#F2B701',
                       'Mac': '#80BA5A', 
                       'NPC': '#7F3C8D', 
                       'Neuron': '#4b4b8f', 
                       'OPC': '#f97b72', 
                       'Oligo': '#ACA4E2', 
                       'Prolif.Metab': '#E73F74', 
                       'Reactive.Ast': '#D4D915',
                       'Vasc': '#CF1C90'}

# hallmark related variables
hallmarks = {
  'HALLMARK_ADIPOGENESIS': 'Adipogenesis', 
  'HALLMARK_ALLOGRAFT_REJECTION': 'Allograft Rejection',
  'HALLMARK_ANDROGEN_RESPONSE': 'Androgen Response', 
  'HALLMARK_ANGIOGENESIS': 'Angiogenesis', 
  'HALLMARK_APICAL_JUNCTION': 'Apical Junction', 
  'HALLMARK_APICAL_SURFACE': 'Apical Surface', 
  'HALLMARK_APOPTOSIS': 'Apoptosis', 
  'HALLMARK_BILE_ACID_METABOLISM': 'Bile Acid Metabolism',
  'HALLMARK_CHOLESTEROL_HOMEOSTASIS': 'Cholesterol Homeostasis', 
  'HALLMARK_COAGULATION': 'Coagulation',
  'HALLMARK_COMPLEMENT': 'Complement',
  'HALLMARK_DNA_REPAIR': 'DNA Repair', 
  'HALLMARK_E2F_TARGETS': 'E2F Targets', 
  'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION': 'EMT', 
  'HALLMARK_ESTROGEN_RESPONSE_EARLY': 'Estrogen Response Early', 
  'HALLMARK_ESTROGEN_RESPONSE_LATE': 'Estrogen Response Late', 
  'HALLMARK_FATTY_ACID_METABOLISM': 'Fatty Acid Metabolism', 
  'HALLMARK_G2M_CHECKPOINT': 'G2M Checkpoint',
  'HALLMARK_GLYCOLYSIS': 'Glycolysis', 
  'HALLMARK_HEDGEHOG_SIGNALING': 'Hedgehog Signaling',
  'HALLMARK_HEME_METABOLISM': 'Heme Metabolism', 
  'HALLMARK_HYPOXIA': 'Hypoxia',
  'HALLMARK_IL2_STAT5_SIGNALING': 'IL2-STAT5 Signaling',
  'HALLMARK_IL6_JAK_STAT3_SIGNALING': 'IL6/JAK/STAT3 Signaling',
  'HALLMARK_INFLAMMATORY_RESPONSE': 'Inflammatory Response', 
  'HALLMARK_INTERFERON_ALPHA_RESPONSE': 'INFA Response', 
  'HALLMARK_INTERFERON_GAMMA_RESPONSE': 'INFG Response',
  'HALLMARK_KRAS_SIGNALING_DN': 'KRAS Signaling Down',
  'HALLMARK_KRAS_SIGNALING_UP': 'KRAS Signaling Up', 
  'HALLMARK_MITOTIC_SPINDLE': 'Mitotic Spindle', 
  'HALLMARK_MTORC1_SIGNALING': 'mTORC1 Signaling',
  'HALLMARK_MYC_TARGETS_V1': 'MYC Targets v1',
  'HALLMARK_MYC_TARGETS_V2': 'MYC Targets v2',
  'HALLMARK_MYOGENESIS': 'Myogenesis', 
  'HALLMARK_NOTCH_SIGNALING': 'Notch Signaling', 
  'HALLMARK_OXIDATIVE_PHOSPHORYLATION': 'Oxidative phosphorylation', 
  'HALLMARK_P53_PATHWAY': 'p53', 
  'HALLMARK_PANCREAS_BETA_CELLS': 'Pancreas Beta Cells', 
  'HALLMARK_PEROXISOME': 'Peroxisome', 
  'HALLMARK_PI3K_AKT_MTOR_SIGNALING': 'PI3K/AKT/mTOR Signaling', 
  'HALLMARK_PROTEIN_SECRETION': 'Protein Secretion', 
  'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY': 'ROS',  
  'HALLMARK_SPERMATOGENESIS': 'Spermatogenesis',
  'HALLMARK_TGF_BETA_SIGNALING': 'TGFB Signaling',
  'HALLMARK_TNFA_SIGNALING_VIA_NFKB': 'TNFA Signaling via NFKB',
  'HALLMARK_UNFOLDED_PROTEIN_RESPONSE': 'UPR',
  'HALLMARK_UV_RESPONSE_DN': 'UV Response Down',
  'HALLMARK_UV_RESPONSE_UP': 'UV Response Up', 
  'HALLMARK_WNT_BETA_CATENIN_SIGNALING': 'Wnt/beta-catenin Signaling', 
  'HALLMARK_XENOBIOTIC_METABOLISM': 'Xenobiotic Metabolism'}

markers_1 = [
    "SNAP25", "SYT1", "SLC17A7", # Neuron
    "VWF", "PDGFRB", "ESAM", # Vasc
    "VEGFA", "NDRG1", "PGK1", # MES-Hyp (malig.)
    "CD14", "CD74", "FCGR3A", # Mac
    "S100B", "AQP1", "TTYH1", # AC (malig.)
    "SOX10", "BCAN", "BCAS1", # OPC (malig.)
    "PLP1", "MAG", "CLDND1", # Oligo
    "ATFX", "NFX1", "CCNI", "KMT2E", "SOX4", # Chromatin-Reg (malig.)
    # "MALAT1", "MT-CO3", "MT-ATP8", # low quality (excluded)
    "EFEMP1", "CHI3L1", "CD44", # MES (malig.)
    "CCND1", "SLC16A1", "CTNNB1", "SOX9", # Prolif-Metab (malig.)
    "GAP43", "KCNF1", "HOXA9", "FABP7", # MES-Ast (malig.)
    "AQP4", "GJA1", "GSN", "HLA-DRA", # Reactive-Ast
    "INA", "NELL2", "DCX", "STMN2", # NPC (malig.)
    "CCL2", "IL1B", "CCL3", "STAT1" # Inflammatory-Mac
] 

markers_2 = [
    "GFAP", "AQP4", "ALDH1L1", "SLC1A3", "CLU", "S100B",
    "PDGFRA", "OLIG1", "OLIG2", "SOX10", "CSPG4",
    "SOX2", "NES", "ASCL1", "DLL3",
    "CHI3L1", "CD44", "ANXA1", "LGALS3", "FN1", "VIM",
    "CA9", "CA12", "VEGFA", "LDHA", "ENO1", "ADM",
    "MKI67", "TOP2A", "CDK1", "HMGB2",
    "P2RY12", "CX3CR1", "TMEM119", "AIF1",
    "CD68", "C1QC", "CTSD", "LYZ", "LGALS3",
    "CD3D", "CD3E",
    "CD8A", "CD8B", "GZMB", "PRF1",
    "FOXP3", "IL2RA",
    "NKG7", "KLRD1", "GNLY",
    "PECAM1", "VWF", "KDR", "ESAM",
    "RGS5", "PDGFRB", "CSPG4", "MCAM",
    "COL1A1", "COL3A1", "DCN", "LUM"
]

marker_genes = list(
    set(markers_1) | 
    set(markers_2) )

def split_at_natural_break(s, max_len=40):
    if len(s) <= max_len:
        return s
    
    # Try to split at the last occurrence of a delimiter before max_len
    delimiters = ['|', '_', '-', ' ']
    for delimiter in delimiters:
        last_pos = s.rfind(delimiter, 0, max_len + 1)  # Find last delimiter before max_len
        if last_pos > 0:  # If found, split there
            return s[:last_pos] + '\n' + s[last_pos + 1:]
    
    # If no delimiter, split at max_len (hard break)
    return s[:max_len] + '\n' + s[max_len:]

def filter_genes(adata, min_count = 1, min_pct = 0.05):
    X = adata.raw.X if adata.raw is not None else adata.X
    if not isinstance(X, np.ndarray):
        X = X.toarray()

    min_cells = int(min_pct * adata.n_obs)
    gene_mask = (X > min_count).sum(axis=0) >= min_cells
    selected_genes = adata.var_names[gene_mask]
    return sorted(selected_genes)