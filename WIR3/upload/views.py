from django.shortcuts import render
import os
from django.http import JsonResponse, HttpResponseForbidden, HttpResponseBadRequest
#import scanpy as sc
#import anndata
#from scipy.sparse import csr_matrix


# Create your views here.

def render_base(request):
    return render(request, "upload/base.html")


def upload(pathname):

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
