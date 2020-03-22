from django.shortcuts import render
import scanpy as sc
import anndata
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.io import loadmat
from scipy.sparse import csr_matrix

# Create your views here.
# Assume all the input here is a scanpy anndata object

def qc_SciPipe(data):
    print("Performing  scipipe for quality control...")
    import seaborn as sns
    sc.pp.calculate_qc_metrics(data, inplace=True)

def normalization_Linnorm(adata):
    print("Performing Linnorm for normalization...")
    normalize(adata)

def normalization_Scran(adata):
    print("Performing Scran for normalization...")
    normalize(adata)

def normalization_Scone(adata):
    print("Performing Scone for normalization...")
    normalize(adata)

def dimReduction_PCA(adata):
    print("Performing PCA for dimension reduction...")
    dimension_reduction(adata)

def dimReduction_tSNE(adata):
    print("Performing tSNE for dimension reduction...")
    dimension_reduction(adata)

def normalize(adata):
    normalized_adata = sc.pp.normalize_total(adata, target_sum=1)
    return normalized_adata

def dimension_reduction(data):
    #get true labels
    m,n = data.obs.shape
    if n >0:
        labels = pd.unique( data.obs.iloc[:,n-1])
    else:
        labels = pd.unique( data.obs.index)
    if len(labels)!=0:
        num_cluster = len(labels)
    else:
        num_cluster = m

    #TSNE
    from openTSNE import TSNE
    tsne_embedded = TSNE().fit(data.X)

    #UMAP
    import umap
    umap_embedded = umap.UMAP(n_neighbors=5,
                          min_dist=0.3,
                          metric='correlation').fit_transform(data.X)

    return tsne_embedded, umap_embedded
