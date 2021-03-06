from django.shortcuts import render
import os
from django.http import JsonResponse, HttpResponseForbidden, HttpResponseBadRequest
from django.views.decorators.csrf import csrf_exempt
from django.core.files.storage import FileSystemStorage
from django.conf import settings
import preprocess.views as preprocess
import analysis.views as analysis
import anndata

# Create your views here.

@csrf_exempt
def render_base(request):
    if request.method == 'POST':
        #checking what the request is for
        if 'qcOption' in request.POST:
            readPath = os.path.join(settings.BASE_DIR, settings.MEDIA_ROOT, 'adata.h5ad')
            data = anndata.read_h5ad(readPath)
            options = request.POST.items()
            for key, value in options:
######### Call Preprocess & Analysis Functions
                if value == 'ScPipe':
                    preprocess.qc_SciPipe(data)
                if value == 'Linnorm':
                    if key == 'normOption':
                        preprocess.normalization_Linnorm(data)
                    else:
                        analysis.diffExpression_Linnorm(data)
                if value == 'Scran':
                    preprocess.normalization_Scran(data)
                if value == 'Scone':
                    preprocess.normalization_Scone(data)
                if value == 'PCA':
                    if key == 'dimReducOption':
                        preprocess.dimReduction_PCA(data)
                    else:
                        analysis.trajectInference_PCA()
                if value == 'tSNE':
                    preprocess.dimReduction_tSNE(data)
                if value == 'Seurat_0.6':
                    analysis.clustering_seurat_06(data)
                if value == 'Seurat_pipe':
                    analysis.clustering_seurat_Pipe(data)
                if value == 'RacelD3':
                    analysis.clustering_seurat_RaceID3(data)
                if value == 'KNN-smoothing':
                    analysis.diffExpression_KNNSmoothing(data)
                if value == 'Linnorm+Drimpute':
                    analysis.diffExpression_LinnormDrimpute(data)
                if value == 'scMerge_s':
                    analysis.cellAnnotation_scMerge_s()
                if value == 'seMerge_us':
                    analysis.cellAnnotation_seMerge_us()
                if value == 'MNNs':
                    analysis.cellAnnotation_MNNs()
                if value == 'Monocle2_DDRTree':
                    analysis.trajectInference_Monocle_DDRTree()
                if value == 'TSCAN':
                    analysis.trajectInference_Tscan()
        else:
            uploaded_file = request.FILES['doc']
            fs = FileSystemStorage()
            fs.save(uploaded_file.name, uploaded_file)
            pathname = os.path.join(settings.BASE_DIR, settings.MEDIA_ROOT, uploaded_file.name)

            print(pathname)
            data = upload(pathname)

            savePath = os.path.join(settings.BASE_DIR, settings.MEDIA_ROOT, 'adata.h5ad')

            print(data.shape)
            data.write(savePath)
    return render(request, "upload/base.html")

def download(email):
    print("link sent to ", email)
    #TODO
    return

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
