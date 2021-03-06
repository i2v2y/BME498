from django.shortcuts import render
import sys
import scanpy as sc
import scipy
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Create your views here.
# Assume all the input here is a scanpy anndata object

def clustering_seurat_06(adata):
    print("Performing seurat_06 for clustering...")
    clustering(adata)

def clustering_seurat_Pipe(adata):
    print("Performing seurat_Pipe for clustering...")
    clustering(adata)

def clustering_seurat_RaceID3(adata):
    print("Performing RaceID3 for clustering...")
    clustering(adata)

def diffExpression_Linnorm(data):
    print("Performing Linnorm for differential gene analysis...")
    sdata = DE(data,10, gene_labels = None)

def diffExpression_KNNSmoothing(data):
    print("Performing KNN Smoothing for differential gene analysis...")
    sdata = DE(data,10, gene_labels = None)

def diffExpression_LinnormDrimpute(data):
    print("Performing Linnorm + Drimpute for differential gene analysis...")
    sdata = DE(data,10, gene_labels = None)

def cellAnnotation_scMerge_s():
    print("Performing scMerge_s for cell annotation...")

def cellAnnotation_seMerge_us():
    print("Performing seMerge_us for cell annotation...")

def cellAnnotation_MNNs():
    print("Performing MNNs for cell annotation...")

def trajectInference_PCA():
	print("Performing PCA for trajectory inference...")

def trajectInference_Monocle_DDRTree():
	print("Performing Monocle DDRTree for trajectory inference...")

def trajectInference_Tscan():
	print("Performing Tscan for trajectory inference...")


def clustering(adata):
    print("Analyzation: Executing clustering...")
    sc.settings.set_figure_params(dpi=80)
    sc.pp.normalize_total(adata, target_sum=1e4)
    adata.X = scipy.sparse.csr_matrix(adata.X)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    return adata


def DE(dataframe, topn, gene_labels = None):

    data = dataframe.X.transpose()

    #get true labels
    m,n = dataframe.obs.shape
    if n >0:
        labels = dataframe.obs.iloc[:,n-1].to_numpy()
        uni_labels = pd.unique( dataframe.obs.iloc[:,n-1])
    else:
        labels = np.ones(m)
        uni_labels = pd.unique( labels)
        topn = 50
    num_cluster = len(uni_labels)


    print("Analyzation: Executing Differential Gene Analysis...")
    #Identification of marker genes and plots gene-cell heatmap of topn markers
    from scipy.io import loadmat
    from scipy.sparse import csr_matrix
    import numpy as np
    from sklearn.cluster import KMeans
    import importlib
    import matplotlib.pyplot as plt
    #TODO: add labels for the genes
    #data gene expression matrix gene x cell
    num_gene, num_cell = data.shape
    cluster_order = []
    gene_idxv = []
    num_cells_inC = []

    #calculate mean of gene expression for each cluster
    #gene x cluster
    gene_mean = np.zeros((num_gene,num_cluster))
    gene_DE_score = np.zeros((num_gene,num_cluster))
    for i in range(num_cluster):
        ind = np.where(labels.flatten() == np.unique(labels)[i])
        gene_mean[:,i] = np.mean(data[:,ind[0]], axis=1).flatten()
        cluster_order = np.append(cluster_order,ind).astype(int)
        num_cells_inC = np.append(num_cells_inC, len(ind[0])).astype(int)
    #find out which cluster expressed the most for each gene
    gene_value_idx = np.argmax(gene_mean,axis = 1)
    #compute DE score for each gene
    for i in range(num_cluster):
        diff = abs( np.matmul(gene_mean[:,i].reshape((num_gene,1)),np.ones((1,num_cluster))) - gene_mean)
        gene_DE_score[:,i] = np.sum( diff, axis = 1)

    #top k for each cluster based on DE score
    gclusters = []
    gscore = []
    for i in range(num_cluster):
        zz_idx = np.where(gene_value_idx == i ) # find genes' DE score in cluster i
        zz_DEscore = gene_DE_score[zz_idx,i].flatten() * (-1)
        zzvalue = np.sort(zz_DEscore, axis = 0) * (-1)
        zz1 = np.argsort(zz_DEscore, axis = 0).astype(int)
        gene_idxv = np.append(gene_idxv, zz_idx[0][zz1[0:topn]] ).astype(int)
        gclusters = np.append(gclusters, i*np.ones((topn,1)) )
        gscore = np.append(gscore, zzvalue[0:topn] )
    GL500 = np.array([gene_idxv,gclusters,gscore]).transpose()

    #datav number of cluster * top n x number of cells
    #column vector, for each cell, the top 10 genes for each cluster that
    #differentiats the most
    datav = data[gene_idxv,:][:,cluster_order]
    #center is num_cluster * top n column vector, each value corresponds the average expression
    #for a gene across all num_cells
    center = np.mean(datav, axis = 1)
    #standard deviation
    scale = np.std(datav, axis = 1)
    #Check for zeros and set them to 1 so not to scale them.
    scale_ind = np.where(scale == 0)
    scale[scale_ind] = 1
    #(data - mean)/scale deviation
    sdata = np.divide(np.subtract(datav,center.reshape((num_cluster*topn,1))),scale.reshape((num_cluster*topn,1)))
    thresh = 3
    return sdata
