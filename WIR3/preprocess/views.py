from django.shortcuts import render

# Create your views here.
# Assume all the input here is a scanpy anndata object

def qc_SciPipe():
    print("Performing  scipipe for quality control...")

def normalization_Linnorm():
    print("Performing Linnorm for normalization...")

def normalization_Scran():
    print("Performing Scran for normalization...")

def normalization_Scone():
    print("Performing Scone for normalization...")

def dimReduction_PCA():
    print("Performing PCA for dimension reduction...")

def dimReduction_tSNE():
    print("Performing tSNE for dimension reduction...")
