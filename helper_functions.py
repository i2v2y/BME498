#!/usr/bin/env python
# coding: utf-8

# In[40]:


def upload(pathname):
    import scanpy as sc
    import os
    import anndata
    from scipy.sparse import csr_matrix
    filename, file_extension = os.path.splitext(pathname)
    if file_extension == ".mat":
        from scipy.io import loadmat
        import pandas as pd
        x = loadmat(pathname)
        keys = []
        for key in x.keys():
            keys.append(key)

        #obs is the cell
        #var is gene
        #pick the largest
        largest = 3
        largest_size = 0
        for i in range(len(keys) - 3):
            if len(x[keys[i+3]].shape) == 2 :
                size = ( x[keys[i+3]].shape[0] * x[keys[i+3]].shape[1] )
            else:
                size = x[keys[i+3]].shape[0]
            if size >= largest_size:
                largest = i+3
                largest_size = size
        obs_d,var_d = {},{}
        for i in range(len(keys) - 3):
            if i != largest-3:
                if (x[keys[i+3]].flatten() ).shape[0] == (x[keys[largest]]).shape[0]:
                    obs_d[keys[i+3]] = x[keys[i+3]].flatten()
                elif (x[keys[i+3]].flatten()).shape[0] == (x[keys[largest]]).shape[1]:
                    var_d[keys[i+3]] = x[keys[i+3]].flatten()
                #else:
        obs_df = pd.DataFrame(data=obs_d)
        var_df = pd.DataFrame(data=var_d)
        
        data = anndata.AnnData(X=x[keys[largest]].todense(),obs=None if obs_df.empty else obs_df,var=None if var_df.empty else var_df)
        
    elif file_extension == ".npz":
        import numpy as np 
        import pandas as pd
        x = np.load(pathname)
        #pick largest size file
        largest = 0
        largest_size = 0
        for i in range(len(x.files)):
            if len(x[x.files[i]].shape) == 2 :
                size = ( x[x.files[i]].shape[0] * x[x.files[i]].shape[1] )
            else:
                size = x[x.files[i]].shape[0]
            if size >= largest_size :
                largest = i
                largest_size = size
        obs_d,var_d = {},{}        
        for i in range(len(x.files)):
            if i != largest:
                if len(x[x.files[i]].flatten() ) == len(x[x.files[largest]]):
                    obs_d[x.files[i]] = x[x.files[i]].flatten()
                elif len(x[x.files[i]].flatten() ) == len(x[x.files[largest]][0]):
                    var_d[x.files[i]] = x[x.files[i]].flatten()
                #else:
        obs_df = pd.DataFrame(data=obs_d)
        var_df = pd.DataFrame(data=var_d)
        data = anndata.AnnData(X=x[x.files[largest]],obs=None if obs_df.empty else obs_df,var=None if var_df.empty else var_df)
    elif file_extension == ".mtx":
        data = sc.read_10x_mtx(os.path.dirname(pathname))
        data.X = data.X.todense()
    elif file_extension == ".csv":
        data = sc.read_csv(pathname)
    elif file_extension == ".xlsx":
        data = sc.read_excel(pathname)
    elif file_extension == ".txt":
        data = sc.read_text(pathname)
    else:
        data = sc.read(pathname)
    
    
    print(pathname, " uploaded !")
    return data

def download(email):
    print("link sent to ", email)
    #TODO
    return

def quality_control():
    print("Preprocessing: Executing QC...")
    return

def dimension_reduction(data):
    print("Preprocessing: Executing Dimension Reduction...")
    #get true labels
    _,n = data.obs.shape
    if n > 1:
        true_labs = data.obs.iloc[:,n-1]
    else :
        true_labs = data.obs
    
    import matplotlib.pyplot as plt
    #TSNE
    from openTSNE import TSNE
    tsne_embedded = TSNE().fit(data.X)
    fig = plt.figure( figsize=(16,7) )
    plt.scatter(tsne_embedded[:, 0], tsne_embedded[:, 1], c=true_labs, s=1.5, cmap='Spectral')
    plt.title(('t-SNE visualization'))
    
    #UMAP
    import umap
    umap_embedded = umap.UMAP(n_neighbors=5,
                          min_dist=0.3,
                          metric='correlation').fit_transform(data.X)
    fig = plt.figure( figsize=(16,7) )
    plt.scatter(umap_embedded[:, 0], umap_embedded[:, 1], c=true_labs, s=1.5, cmap='Spectral')
    plt.title('UMAP visualization')
    
    return tsne_embedded, umap_embedded

def alignment():
    print("Preprocessing: Executing alignment...")
    return

def normalization():
    print("Preprocessing: Executing normalization...")
    return

def clustering():
    print("Analyzation: Executing clustering...")
    return

def DE(data, labels, topn, num_cluster, gene_labels = None):
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
    
    #plot color map
    fig = plt.figure( figsize=(16,16) )
    plt.imshow(sdata,cmap = 'RdBu_r' ,vmin= -thresh, vmax=  thresh, aspect='auto')
    #plt.colorbar()
    xtkval = np.cumsum(num_cells_inC)
    xtkval1 = np.zeros(num_cluster)
    xtllab = []
    for i in range(num_cluster):
        #plt.axhline(y=i*topn, color='w')
        if i == 0:
            xtkval1[i] = 0.5*num_cells_inC[i]
            xtllab = np.append(xtllab, "C1")
            plt.axhline(xmin=0, xmax=xtkval[i]/num_cell, y=num_cluster*topn, linewidth=4, color='C'+str(i))
        else:
            xtkval1[i] = 0.5*num_cells_inC[i] +  xtkval[i-1]
            xtllab = np.append(xtllab, "C"+str(i+1))
            plt.axhline(xmin=xtkval[i-1]/num_cell, xmax=xtkval[i]/num_cell, y=num_cluster*topn,linewidth=4, color='C'+str(i%9))
    plt.xticks(xtkval1,xtllab)
    if gene_labels is None:
        plt.yticks(np.arange(num_cluster*topn))
    else:
        plt.yticks(np.arange(num_cluster*topn),gene_labels.flatten())
    #filename="marker_gene_" + str(num_cluster) + "_" +str(topn) + "_" + str(num_cell) + ".eps"
    #plt.savefig(filename)
    plt.show()
    return sdata

def TI():
    print("Analyzation: Executing Trajectory Inference...")
    return

def cell_annotation():
    print("Analyzation: Executing cell_annotation...")
    return

def visualize():
    print("Plotting")
    return

def preprocess(feature, userid, data ):
    if feature == "quality_control":
        quality_control()
    elif feature == "dimension_reduction": 
        dimension_reduction(data)
    elif feature == "alignment": 
        alignment() 
    elif feature == "normalization": 
        normalization() 
    else:
        print("Error, feature %s not supported" %(feature))
        
    visualize()
    
    return

def analyze(feature, userid ):
    if feature == "clustering":
        clustering()
    elif feature == "DE": 
        #get true labels
        _,n = data.obs.shape
        if n > 1:
            true_labs = data.obs.iloc[:,n-1].to_numpy()
        else :
            true_labs = data.obs.to_numpy()
        num_cluster = len(np.unique(true_labs))
        sdata = DE(data.X.transpose(), true_labs, 10, num_cluster, gene_labels = None)
    elif feature == "TI": 
        TI() 
    elif feature == "cell_annotation": 
        cell_annotation() 
    else:
        print("Error, feature %s not supported" %(feature))
        
    visualize()
    return

def login(userid, password):
    print("user %s successfully logged in" %(userid))
    return


# In[2]:


def parse_request(data):
    for i in range (len(data)):
        userid = data[i][0]
        option = data[i][1]
        detail = data[i][2]
        #print(userid, option, detail)
        if option == "upload":
            upload(detail)
        elif option == "download":
            download(detail)
        elif option == "preprocess":
            preprocess(detail,userid)
        elif option == "analyze":
            analyze(detail,userid)
        elif option == "login":
            login(userid, detail)
        else:
            print("Unable to parse request: %s" %(data[i]))
    return
    
    
    


# In[5]:


import csv

with open('test.csv', newline='') as csvfile:
    data = list(csv.reader(csvfile))

#print(data)
parse_request(data)


# In[59]:


#filename = "/Users/yumengluo/Desktop/claire/research/Test_5_Zeisel.mat"
filename = '/Users/yumengluo/Desktop/claire/research/cortex.npz'
#filename = '/Users/yumengluo/Desktop/claire/BME498/data/filtered_gene_bc_matrices/hg19/matrix.mtx'

data = upload(filename)


# In[39]:


dimension_reduction(data)


# In[60]:


#get true labels
_,n = data.obs.shape
if n > 1:
    true_labs = data.obs.iloc[:,n-1].to_numpy()
else :
    true_labs = data.obs.to_numpy()
num_cluster = len(np.unique(true_labs))
DE(data.X.transpose(), true_labs, 10, num_cluster, gene_labels = None)


# In[53]:


data.X.shape


# In[ ]:





# In[ ]:



