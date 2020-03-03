from django.shortcuts import render

# Create your views here.
# Assume all the input here is a scanpy anndata object

def qc_SciPipe():
    print("p_scipipe()")

def alignment_RsuBread():
    print("p_RsuBread()")

def normalization_Linnorm():
    print("p_Linnorm()")

def normalization_Scran():
    print("p_scran()")

def normalization_Scone():
    print("p_Scone()")

def dimReduction_PCA():
    print("p_pca()")

def dimReduction_tSNE():
    print("p_tSNE()")
