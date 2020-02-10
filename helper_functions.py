#!/usr/bin/env python
# coding: utf-8

# In[147]:


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
        
        data = anndata.AnnData(X=x[keys[largest]],obs=None if obs_df.empty else obs_df,var=None if var_df.empty else var_df)
        
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

def dimension_reduction():
    print("Preprocessing: Executing Dimension Reduction...")
    return

def alignment():
    print("Preprocessing: Executing alignment...")
    return

def normalization():
    print("Preprocessing: Executing normalization...")
    return

def clustering():
    print("Analyzation: Executing clustering...")
    return

def DE():
    print("Analyzation: Executing Differential Gene Analysis...")
    return

def TI():
    print("Analyzation: Executing Trajectory Inference...")
    return

def cell_annotation():
    print("Analyzation: Executing cell_annotation...")
    return

def visualize():
    print("Plotting")
    return

def preprocess(feature, userid ):
    if feature == "quality_control":
        quality_control()
    elif feature == "dimension_reduction": 
        dimension_reduction()
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
        DE()
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


# In[40]:


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
    
    
    


# In[43]:


import csv

with open('test.csv', newline='') as csvfile:
    data = list(csv.reader(csvfile))

#print(data)
parse_request(data)


# In[166]:


#filename = "/Users/yumengluo/Desktop/claire/research/Test_5_Zeisel.mat"
#filename = '/Users/yumengluo/Desktop/claire/research/cortex.npz'
filename = '/Users/yumengluo/Desktop/claire/BME498/data/filtered_gene_bc_matrices/hg19/matrix.mtx'

data = upload(filename)
print(data.X)
print(data.obs)
print(data.var)


# In[167]:


data.obs


# In[168]:


data.var


# In[169]:


data.X


# In[ ]:





# In[ ]:





# In[ ]:




