from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq


def get_activity(adata, key='tfa_stan'):
    s = key[:2]
    adata_tfa = AnnData(
        X = adata.obsm[key],
        obs = adata.obs,
        obsm = {name: obj for (name, obj) in adata.obsm.items() if s not in name},
        layers = {name: obj for (name, obj) in adata.obsm.items() if s in name})
    adata_tfa.uns = adata.uns
    return adata_tfa


def compute_spatial_expression(adata, genes):
    sq.gr.spatial_neighbors(adata, n_rings=1)
    A = adata.obsp['spatial_connectivities']
    mat = sc.get.obs_df(adata, genes)
    mat_neighbor = pd.DataFrame(
        (A + np.eye(adata.n_obs)).dot(mat)/(1+A.sum(axis=1)),
        index = mat.index,
        columns = mat.columns)
    return mat_neighbor


def compute_score_by_subtype(adata_list, subtype, act='TF'):
    dfs_tf = dict()
    for sample in adata_list.keys():
        adata = adata_list[sample]
        if subtype in adata.obs.columns:
            categories = adata.obs[subtype].cat.categories
            mean_df = adata.to_df().groupby(adata.obs[subtype]).mean()
            df = []
            for i in range(len(adata.uns['rank_genes_groups']['names'])):
                for j, cat in enumerate(categories):
                    tf = adata.uns['rank_genes_groups']['names'][i][j]
                    pval = adata.uns['rank_genes_groups']['pvals_adj'][i][j]
                    mean = mean_df.loc[cat, tf]
                    df.append([sample, tf, cat, mean, pval, -np.log10(pval+1e-10), i])

            dfs_tf[sample] = pd.DataFrame(df, columns=['sample', act, 'subtype', act+'a', 'p_adj', 'neg_log_p_adj', 'rank'])
            dfs_tf[sample]['abs_{}a'.format(act)] = dfs_tf[sample][act+'a'].abs()

    dfs_tf_diff = dict()
    for sample in dfs_tf.keys():
        adata = adata_list[sample]
        df = dfs_tf[sample]
        df['diff'] = 0
        mat = adata.to_df()
        for i in df.index:
            tf = df[act][i]
            pat = df['subtype'][i]
            x = mat.loc[adata.obs[subtype]==pat, tf].mean(numeric_only=True)
            y = mat.loc[adata.obs[subtype]!=pat, tf].mean(numeric_only=True)
            df['diff'][i] = x-y
        dfs_tf_diff[sample] = df
        
    df_tf_diff = pd.concat(dfs_tf_diff, ignore_index=True)
    return df_tf_diff


def compute_morans_i(adatas_dict):
    all_tfs = set()
    for adata in adatas_dict.values():
        all_tfs.update(adata.var_names)
    moran_df = pd.DataFrame(index=sorted(all_tfs), columns=list(adatas_dict.keys()))
    
    for sample, adata in adatas_dict.items():
        if 'spatial' in adata.obsm:
            try:
                sq.gr.spatial_neighbors(adata, coord_type='generic', key_added='spatial')
                sq.gr.spatial_autocorr(
                    adata,
                    mode='moran',
                    n_perms=100,  
                    n_jobs=1,     
                    genes=adata.var_names.tolist()  
                )
                
                if 'moranI' in adata.uns:
                    moran_results = adata.uns['moranI']
                    moran_values = moran_results['I']
                    valid_genes = moran_values.index.intersection(adata.var_names)
                    moran_df.loc[valid_genes, sample] = moran_values.loc[valid_genes].values
                else:
                    print(f"Warning: 'moranI' not found in adata.uns for sample {sample}")
                    moran_df.loc[adata.var_names, sample] = np.nan
                    
            except Exception as e:
                print(f"Error computing Moran's I for sample {sample}: {e}")
                moran_df.loc[adata.var_names, sample] = np.nan
        else:
            print(f"Warning: No spatial coordinates found for sample {sample}")
            moran_df.loc[adata.var_names, sample] = np.nan
    
    return moran_df


def compute_bivariate_morans_i(adatas_var1, adatas_var2, is_symmetric=False, keep_diagonal=False):
    samples = list(set(adatas_var1.keys()) & set(adatas_var2.keys()))
    if not samples:
        raise ValueError("No common samples found")
    
    all_var1s = sorted(set().union(*[set(adatas_var1[s].var_names) for s in samples]))
    all_var2s = sorted(set().union(*[set(adatas_var2[s].var_names) for s in samples]))
    
    # Generate pairs based on symmetric and keep_diagonal flags
    var1_var2_pairs = []
    for var1 in all_var1s:
        for var2 in all_var2s:
            # Skip diagonal if keep_diagonal is False
            if not keep_diagonal and var1 == var2:
                continue
            # For symmetric mode, only keep var1 <= var2 (alphabetical order)
            if is_symmetric:
                if var1 <= var2:  # This ensures we only get each unordered pair once
                    var1_var2_pairs.append(f"{var1}_{var2}")
            else:
                # In non-symmetric mode, keep all combinations
                var1_var2_pairs.append(f"{var1}_{var2}")
    
    bivariate_df = pd.DataFrame(index=var1_var2_pairs, columns=samples, dtype=float)
    
    for sample in samples:
        adata_var1 = adatas_var1[sample]
        adata_var2 = adatas_var2[sample]
        
        if 'spatial' not in adata_var1.obsm:
            print(f"Warning: No spatial coordinates for {sample}")
            bivariate_df[sample] = np.nan
            continue
        
        try:
            if 'spatial_connectivities' not in adata_var1.obsp:
                sq.gr.spatial_neighbors(adata_var1, key_added='spatial')
            
            var1s_in_sample = list(adata_var1.var_names)
            var2s_in_sample = list(adata_var2.var_names)
            y_var1 = adata_var1[:, var1s_in_sample].X.toarray()  # n_cells × n_var1s
            y_var2 = adata_var2[:, var2s_in_sample].X.toarray()  # n_cells × n_var2s
            
            y_var1_centered = y_var1 - y_var1.mean(axis=0, keepdims=True)
            y_var2_centered = y_var2 - y_var2.mean(axis=0, keepdims=True)
            W = adata_var1.obsp['spatial_connectivities'].toarray()
            W_sum = W.sum()
            
            if W_sum == 0:
                print(f"Warning: Zero weight sum for {sample}")
                bivariate_df[sample] = np.nan
                continue
            
            std_var1 = y_var1_centered.std(axis=0, ddof=1)
            std_var2 = y_var2_centered.std(axis=0, ddof=1)
            std_var1[std_var1 == 0] = 1
            std_var2[std_var2 == 0] = 1
            
            z_var1 = y_var1_centered / std_var1  # n_cells × n_var1s
            z_var2 = y_var2_centered / std_var2  # n_cells × n_var2s
            # I_ij = (z_var1_iᵀ @ W @ z_var2_j) / W_sum
            bivariate_matrix = (z_var1.T @ W @ z_var2) / W_sum
            
            for i, var1 in enumerate(var1s_in_sample):
                for j, var2 in enumerate(var2s_in_sample):
                    # Apply same filtering rules when populating the dataframe
                    if not keep_diagonal and var1 == var2:
                        continue
                    
                    if is_symmetric:
                        if var1 <= var2:
                            pair_name = f"{var1}_{var2}"
                            bivariate_df.loc[pair_name, sample] = bivariate_matrix[i, j]
                    else:
                        pair_name = f"{var1}_{var2}"
                        bivariate_df.loc[pair_name, sample] = bivariate_matrix[i, j]
            
        except Exception as e:
            print(f"Error computing for sample {sample}: {e}")
            bivariate_df[sample] = np.nan
    
    return bivariate_df
