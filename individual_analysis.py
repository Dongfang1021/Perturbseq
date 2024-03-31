#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as plt
import gseapy
import pertpy as pt
import muon as mu
from matplotlib.pyplot import figure

import plotnine
import string
import argparse


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument('--h5ad_input', help='gRNA file path')
argv = vars(parser.parse_args())

if argv['h5ad_input'] == None:
    	raise Exception ('You should provide h5ad file!')
else:
	h5ad_data=argv['h5ad_input'].strip()



# In[16]:

gene_symbol = h5ad_data.split('.h5ad')[0]
data = sc.read_h5ad(h5ad_data)



def basic_stat(data, gene_symbol):
		#mkdir folder for each sample
	gene_list = ['IKBKG',
			'GRHL2',
			'REEP3',
			'TAB3',
			'CHUK',
			'AMOTL2',
			'PPP2CA',
			'TAB2',
			'TNFRSF1A',
			'IRF6',
			'NUP62CL',
			'RAB7A',
			'REEP5',
			'RDH14',
			'RELA',
			'NFKB1',
			'TRADD',
			'RBL2',
			'PTEN',
			'UBC',
			'TNFRSF1B',
			'TRAF2',
			'C1orf43',
			'LATS2',
			'IKBKB',
			'RB1',
			'TRAF1',
			'PDCD10',
			'IMMP1L',
			'STK17B',
			'SMPD1',
			'MAP3K8',
			'ALAS1',
			'KCTD5',
			'DMD',
			'VCP',
			'MBNL3',
			'ACTG1',
			'GPI',
			'NRIP1',
			'SNRPD3',
			'RIPK1',
			'NFKB2',
			'TAB1']
	os.makedirs(gene_symbol, exist_ok=True)
	sc.pl.dotplot(data, gene_list, groupby="Treatment_CRISPR", dendrogram=True, save=gene_symbol+"dotplot_mean_expression.png")
	sc.pl.stacked_violin(data, gene_list, groupby="Treatment_CRISPR", save=gene_symbol+"violin_mean_expression.png")
	sc.tl.dendrogram(data, groupby='Treatment_CRISPR')
	sc.pl.dendrogram(data, groupby='Treatment_CRISPR', save=gene_symbol+"_dendrogram.png")
	#sc.pl.clustermap(data, groupby='Treatment_CRISPR')
	sc.pl.pca(data, color="Treatment_CRISPR", save=gene_symbol+"_pca.png")
		# whole dataset umap
	sc.tl.rank_genes_groups(data, groupby='Treatment_CRISPR', method="wilcoxon")
	sc.pl.rank_genes_groups_dotplot(data, var_names=[gene_symbol],
        values_to_plot="logfoldchanges", cmap='bwr',
        vmin=-4,
        vmax=4,
        colorbar_title='log fold change', save= '_'+gene_symbol+'_foldchanges_gene_level.png')
	# allcell_name_array = data.uns['rank_genes_groups']['names']
	# data.uns['rank_genes_groups']
	# allcell_pvals_adj_array = data.uns['rank_genes_groups']['pvals_adj']
	# allcell_pval_array = data.uns['rank_genes_groups']['pvals']
	# allcell_log_array = data.uns['rank_genes_groups']['logfoldchanges']
	# allcell_gene_list = allcell_name_array[gene_symbol]
	# allcell_pvals_adj_array = allcell_pvals_adj_array[gene_symbol]
	# allcell_pvals_array = allcell_pval_array[gene_symbol]
	# allcell_logchanges_array = allcell_log_array[gene_symbol]
	# allcell_df = pd.DataFrame(list(zip(allcell_gene_list, allcell_pvals_adj_array, allcell_pvals_array, allcell_logchanges_array)), columns=['gene_symbol', 'pvals_adj','pvals', 'logfoldchange'])
	# allcell_df[allcell_df['gene_symbol']==gene_symbol].to_csv(gene_symbol+"target_DEG.xls", sep='\t')
	# allcell_df_selected = allcell_df[allcell_df['pvals']<0.1]
	# allcell_df_selected.to_csv(gene_symbol+"cells_DEG.xls", sep='\t')
	# top_10_up = list(allcell_df_selected[allcell_df_selected['logfoldchange']>0].sort_values(by=['pvals', 'logfoldchange'], ascending=[True, False])['gene_symbol'])[:10]
	# top_10_down = list(allcell_df_selected[allcell_df_selected['logfoldchange']<0].sort_values(by=['pvals', 'logfoldchange'], ascending=[True, True])['gene_symbol'])[:10]
	# top_genes = top_10_up + top_10_down
	sc.pl.heatmap(data, var_names=gene_list, groupby='Treatment_CRISPR', use_raw=False, cmap='inferno', figsize=(12, 8), dendrogram=True)
	plt.savefig(gene_symbol+'/'+gene_symbol+"_top10_genes_heatmap.png")
	
	sc.pp.normalize_total(data, target_sum=1e4)
	sc.pp.log1p(data)
	sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5)
	data = data[:, data.var.highly_variable]
	sc.pp.regress_out(data, ["total_counts", "pct_counts_mt"])
	sc.pp.scale(data, max_value=10)
	sc.tl.pca(data, svd_solver="arpack")
	sc.pp.neighbors(data,  n_neighbors=30)
	sc.tl.leiden(data, resolution=0.5)
	sc.tl.umap(data)
	sc.pl.umap(data, color=["leiden", "Treatment_CRISPR"], save="_"+gene_symbol+".png")
	sc.tl.tsne(data)
	sc.pl.tsne(data, color=["leiden", "Treatment_CRISPR"], save="_"+gene_symbol+".png")

basic_stat(data, gene_symbol)