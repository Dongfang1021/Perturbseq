#!/usr/bin/env python
# coding: utf-8

# In[13]:


import scanpy as sc
import pandas as pd
import argparse
import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import numpy as np
import seaborn as sns
from anndata.experimental.multi_files import AnnCollection
import anndata as ad 

# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument('--sample', help='sample name list')
parser.add_argument('--feature_library', help='feature library file')
argv = vars(parser.parse_args())

if argv['sample'] == None:
    raise Exception ('You should provide sample name list!')
else:
    sample=argv['sample'].strip()
    
if argv['feature_library'] == None:
    raise Exception ('you should provide feature library file!')
else:
    feature_library = argv['feature_library'].strip()
    
reference = pd.read_csv(feature_library)
ref_dict = dict(zip(list(reference['name']), list(reference['target_gene_name'])))
    


# In[ ]:


# Define the base directory where the subfolders are located
base_dir = '/run/media/Dongfang/DataUSB/Perturb_seq/10XGenomixs_10/'
titles = sample.split(',')
titles = sorted(titles)
h5ad_list = []
for each in titles:
    adata = sc.read_h5ad(each+'.h5ad')
    adata.var_names_make_unique()
    print(adata.obs)
    print(adata.var)
    h5ad_list.append(adata)

adata_merged = h5ad_list[0].concatenate(h5ad_list[1:], batch_key='sample', index_unique=None)
adata_merged.write("merged.h5ad")