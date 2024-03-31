#!/usr/bin/env python
# coding: utf-8

# In[63]:


import numpy as np
import pandas as pd
import scanpy as sc
import os
import matplotlib as plt
import gseapy
import pertpy as pt
import muon as mu
from matplotlib.pyplot import figure
import glob
import anndata as ad
from anndata.experimental.multi_files import AnnCollection


# In[110]:


adata = sc.read_h5ad('merged.h5ad')


# In[112]:


adata = adata[adata.obs['gRNA_count']==1]


# In[113]:


adata.obs


# In[114]:


sc.pl.highest_expr_genes(adata, n_top=20,)


# In[ ]:





# In[115]:


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'


# In[116]:


sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[117]:


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)


# In[118]:


adata = adata[adata.obs.pct_counts_mt < 10, :]


# In[119]:


sc.pl.highest_expr_genes(adata, n_top=50, )


# In[120]:


sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


# In[121]:


sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)


# In[122]:


sc.tl.umap(adata)


# In[123]:


sc.tl.leiden(adata)


# In[124]:


sc.pl.umap(adata, color=['leiden', 'Treatment'])


# In[125]:


adata.obs['NT'] = adata.obs['target_gene'].apply(lambda x: "NT" if x=='NegCtrl' else "Perturbed")


# In[129]:


sc.pl.umap(adata, color=['leiden', 'Treatment', "CHUK", 'IKBKG', 'IRF6', 'RELA', 'TNFRSF1A', 'ALAS1'])


# In[ ]:


adata.obs['perturbation'] = adata.obs['NT'].apply(lambda x: x if x=="NT" else "Perturbed")
adata.obs['NT'] = adata.obs['target_gene']

# In[132]:


adata.write('mix_summary.h5ad')


# In[ ]:

def split_chunks(adata):
    """dict = dict(gene_id_gene_name:guide + NT list)"""
    for each in list(adata.obs['target_gene']):
        if each != "NegCtrl":
              adata_sub = adata[adata.obs['target_gene'].isin([each, "NegCtrl"])]
              adata_sub.write("chunks/"+each+".h5ad")
              
split_chunks(adata)
