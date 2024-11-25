import scanpy as sc
import pandas as pd
import numpy as np
from scipy.spatial.distance import cosine
from sklearn.metrics.pairwise import cosine_similarity as cosine_similarity_sklearn



# cosine similarity
def cosine_similarity(adata, gene_field, embedding):
    embeddings = []
    genes = []

    gene_ids = adata.obs[gene_field].values
    for gene in adata.obs[gene_field].unique():
        genes.append(gene)
        embeddings.append(np.mean(adata.obsm[embedding][gene_ids==gene,:], axis=0))
    
    embedding = np.stack(embeddings)
    genes = np.array(genes)
    return pd.DataFrame(data=cosine_similarity_sklearn(embedding), index=genes, columns=genes)


def centered_embedding_adata(
    adata,
    embeddings_key,
    gene_field,
    guide_field,
    ntc_key,
    pool_ntcs=False
):
    
    """
    Centers an embedding at NTC cells.

    adata: AnnData object
    embeddings_key: key in adata.obsm that contains embedding to be centered, e.g. 'X_pca', 'X_scVI',... Assumed to be a dense
        numpy array. 
    guide_field: field name in adata.obs that contains guide assignment
    gene_field: field name in adata.obs that contains gene-level assignment
    ntc_key: key that corresponds to negative controls in adata.obs[gene_field]
    pool_ntcs: whether to pool NTCs to provide the mean estimate. If False, will first calculate
        means over individual NTC guides and then calculate the overall mean of those.

    Returns the an AnnData object with added centered embeddings_key to adata.obsm.
    """
 
    
    centered_key = embeddings_key+"_centered"
    ntc_idx = adata.obs[gene_field] == ntc_key
    ntc_guides = adata.obs[ntc_idx][guide_field].unique().tolist()
   
    mat = adata.obsm[embeddings_key].copy()
    
    if pool_ntcs:
        mat -= np.mean(mat[ntc_idx], axis=0, keepdims=True)
        
    else:
        means=[]
        for ntc_guide in ntc_guides:
            means.append(np.mean(mat[adata.obs[guide_field].isin(ntc_guides)], axis=0))
        mat -=  np.mean(means, axis=0, keepdims=True)
        
    adata.obsm[centered_key] = mat