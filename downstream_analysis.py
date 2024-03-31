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
# Ensure the 'Report' directory exists
CRISPR_dir = os.path.join(base_dir, 'CRISPR')
os.makedirs(CRISPR_dir, exist_ok=True)

# Copy files into the 'Report' directory and generate a new list with updated paths
individual_files = []
for title in titles:
    original_path = os.path.join(base_dir, title, 'outs', 'crispr_analysis', 'protospacer_calls_per_cell.csv')
    new_path = os.path.join(CRISPR_dir, f"{title}_protospacer_calls_per_cell.csv")  # New file name based on title
    shutil.copy(original_path, new_path)
    individual_files.append({'title': title, 'path': new_path})

QC_dir = os.path.join(base_dir, 'QC')
os.makedirs(QC_dir, exist_ok=True)


QC_files = []
for title in titles:
    original_path = os.path.join(base_dir, title, 'outs', 'metrics_summary.csv')
    new_path = os.path.join(QC_dir, f"{title}_metrics_summary.csv")  # New file name based on title
    shutil.copy(original_path, new_path)
    QC_files.append({'title': title, 'path': new_path})
    
name_list = []
file_list = []
for each in QC_files:
    name_list.append(each['title'])
    file_list.append(pd.read_csv(each['path']))
    QC_sum = pd.concat(file_list)
QC_sum['Treatment'] = name_list
QC_sum.to_csv("QC/QC_summary.xls",sep='\t')



summary_files = []
for title in titles:
    original_path = os.path.join(base_dir, title, 'outs', 'crispr_analysis', 'protospacer_calls_summary.csv')
    new_path = os.path.join(CRISPR_dir, f"{title}_protospacer_calls_summary.csv")  # New file name based on title
    shutil.copy(original_path, new_path)
    summary_files.append({'title': title, 'path': new_path})
def basic_figure_summary(file, name):    
	# Load the CSV file into a Pandas DataFrame
	df = pd.read_csv(file)
	df_total = df.head(4)
	# Assuming 'column1' and 'column2' are numerical columns you want to plot in a scatter plot
	# Creating a bar chart
	colors = ['red', 'blue', 'green', 'orange', 'purple', 'cyan', 'magenta', 'yellow', 'black', 'pink']

	fig, ax = plt.subplots(figsize=(8, 6))  # Adjust the figure size as needed
	bars = ax.bar(df_total['feature_call'], df_total['pct_cells'], color=colors)
	

	# Adding numbers on top of each bar
	for bar in bars:
		height = bar.get_height()
		ax.annotate(f'{height:.2f}',
					xy=(bar.get_x() + bar.get_width() / 2, height),
					xytext=(0, 3),  # 3 points vertical offset
					textcoords="offset points",
					ha='center', va='bottom')
	# Rotate x-axis labels
	ax.set_xticks(range(len(df_total['feature_call'])))
	ax.set_xticklabels(df_total['feature_call'], rotation=45)
	plt.title(name + ' feature call and percentage')
	plt.xlabel('Feature call')
	plt.ylabel('Percentage %')
	plt.tight_layout()
	plt.savefig(base_dir+"CRISPR/"+name+'_feature_call.png')
       

def basic_figure_top10(file, name):    
	# Load the CSV file into a Pandas DataFrame
	df = pd.read_csv(file)
	df_single = df.iloc[4:].sort_values(by='pct_cells', ascending=False).head(10)
	colors = ['red', 'blue', 'green', 'orange', 'purple', 'cyan', 'magenta', 'yellow', 'black', 'pink']
	fig, ax = plt.subplots(figsize=(8, 6))  # Adjust the figure size as needed
	bars = ax.bar(df_single['feature_call'], df_single['pct_cells'], color=colors)
	# Adding numbers on top of each bar
	for bar in bars:
		height = bar.get_height()
		ax.annotate(f'{height:.2f}',
					xy=(bar.get_x() + bar.get_width() / 2, height),
					xytext=(0, 3),  # 3 points vertical offset
					textcoords="offset points",
					ha='center', va='bottom')
	ax.set_xticks(range(len(df_single['feature_call'])))
	ax.set_xticklabels(df_single['feature_call'], rotation=45)
	plt.title(name + ' top10 guides and percentage')
	plt.xlabel('top10 guides')
	plt.ylabel('Percentage %')
	plt.tight_layout()
	plt.savefig(base_dir + 'CRISPR/'+name+'top10_feature_call.png')

	

for each in summary_files:
	basic_figure_summary(each['path'], each['title'])
	basic_figure_top10(each['path'], each['title'])

def basic_figure_total(file, name):    
	# Load the CSV file into a Pandas DataFrame
	df = pd.read_csv(file)
	# Assuming 'column1' and 'column2' are numerical columns you want to plot in a scatter plot
	# Creating a bar chart
	colors = ['red', 'blue', 'green', 'orange', 'purple', 'cyan', 'magenta', 'yellow', 'black', 'pink']

	fig, ax = plt.subplots(figsize=(8, 6))  # Adjust the figure size as needed
	frequency = df['num_features'].value_counts(normalize=True)*100
	frequency_df = frequency.reset_index()
	frequency_df.columns = ['sgRNA per cell', 'Percentage']
	bars = ax.bar(range(len(frequency_df['sgRNA per cell'])), frequency_df['Percentage'], color=colors)
	

	# Adding numbers on top of each bar
	for index, bar in enumerate(bars):
		height = bar.get_height()
		ax.annotate(f'{height:.2f}',
					xy=(index, height),
					xytext=(0, 3),  # 3 points vertical offset
					textcoords="offset points",
					ha='center', va='bottom')
	# Rotate x-axis labels
	ax.set_xticks(range(len(frequency_df['sgRNA per cell'])))
	ax.set_xticklabels(frequency_df['sgRNA per cell'], rotation=45)
	plt.title(name + ' sgRNA number distribution')
	plt.xlabel('sgRNA per cell')
	plt.ylabel('Percentage %')
	plt.tight_layout()
	plt.savefig(base_dir + 'CRISPR/'+name+'_sgRNA_number_distribution.png')
       
def mean_umi(list_umi):
	list_umi = [int(x) for x in list_umi]
	return sum(list_umi)/len(list_umi)

def detected_guide(file, name, reference):    
	# Load the CSV file into a Pandas DataFrame
	df = pd.read_csv(file)
	guide_set = set(reference['name'])
	detected_single_guide_set = set(df[df['num_features']==1]['feature_call'])
	detected_single_guide_umi = df[df['num_features']==1]
	detected_single_guide_umi['num_umis'] = detected_single_guide_umi['num_umis'].apply(int)

	detected_single_guide_umi_mean = df[df['num_features']==1].groupby('feature_call')['num_umis'].apply(list).reset_index(name='mean umi')
	detected_single_guide_umi_mean['mean umi'] = detected_single_guide_umi_mean['mean umi'].apply(mean_umi)

	plt.figure(figsize=(24, 6))
	sns.stripplot(x='feature_call', y='num_umis', data=detected_single_guide_umi, jitter=True, palette='bright', size=3)
	plt.gca().invert_yaxis()
	plt.yscale('log')

	# Save the figure with bbox_inches='tight' to include all labels
	plt.xticks(rotation=90, fontsize=8)
	plt.title(name+'UMI Distribution.png')
	plt.savefig(base_dir + "CRISPR/"+name+'_umi_distribution.png', bbox_inches='tight')
	plt.cla()

	plt.figure(figsize=(24, 6))
	sns.boxplot(x='feature_call', y='num_umis', data=detected_single_guide_umi, palette='bright')
	# Save the figure with bbox_inches='tight' to include all labels
	plt.xticks(rotation=90, fontsize=8)
	plt.yscale('log')
	plt.title(name+'UMI Distribution.png')
	plt.savefig(base_dir + "CRISPR/"+name+'_umi_distribution_boxplot.png', bbox_inches='tight')
	plt.cla()


	plt.figure(figsize=(24, 6))
	sns.lineplot(x='feature_call', y='mean umi', data=detected_single_guide_umi_mean,  palette='bright', size=3)
	plt.yscale('log')
	plt.gca().invert_yaxis()

	# Save the figure with bbox_inches='tight' to include all labels
	plt.xticks(rotation=90, fontsize=8)
	plt.title(name+'Mean UMI Distribution.png')
	plt.savefig(base_dir + "CRISPR/"+name+'_mean_umi_distribution.png', bbox_inches='tight')
	return len(guide_set - detected_single_guide_set)

missing_guide_dict = {}

for each in individual_files:
	basic_figure_total(each['path'], each['title'])
	missing_guide_dict[each['title']] = detected_guide(each['path'], each['title'], reference)

# Creating the bar chart
# The keys of the dictionary will be the categories
categories = list(missing_guide_dict.keys())
# The values of the dictionary will be the heights of the bars
values = list(missing_guide_dict.values())

fig, ax = plt.subplots(figsize=(8, 6))  # Adjust the figure size as needed
colors = ['red', 'blue', 'green', 'orange', 'purple', 'cyan', 'magenta', 'yellow', 'black', 'pink']

bars = ax.bar(range(len(categories)), values, color=colors)
	

	# Adding numbers on top of each bar
for index, bar in enumerate(bars):
	height = bar.get_height()
	ax.annotate(f'{height}',
					xy=(index, height),
					xytext=(0, 3),  # 3 points vertical offset
					textcoords="offset points",
					ha='center', va='bottom')
	# Rotate x-axis labels
ax.set_xticks([x - 1 for x in range(len(categories))])
ax.set_xticklabels(categories, rotation=45)
plt.title('Missing guide barchart')
plt.xlabel('Treatment')
plt.ylabel('Missing guides number')
plt.tight_layout()
plt.savefig(base_dir+"CRISPR/"+'Missing_guide_number.png')


# In[53]:


def process_data(Gene_h5, guide_file, ref_dict, name):
    # Load the expression data and gRNA info
    adata = sc.read_10x_h5(Gene_h5, gex_only=True)
    adata.var_names_make_unique()
    grna_df = pd.read_csv(guide_file)
    barcode_dict = dict(zip(list(grna_df['cell_barcode']), list(grna_df['feature_call'])))
    sgRNA_number_dict = dict(zip(list(grna_df['cell_barcode']), list(grna_df['num_features'])))
    sgRNA_umi_dict = dict(zip(list(grna_df['cell_barcode']), list(grna_df['num_umis'])))
    # Annotate the adata object with the gRNA information
    adata.obs['gRNA'] = [barcode_dict.get(x) for x in list(adata.obs_names)]
    adata.obs['target_gene'] = adata.obs['gRNA'].apply(lambda x: ref_dict.get(x))
    adata.obs['gRNA_umi'] =[sgRNA_umi_dict.get(x) for x in list(adata.obs_names)]
    adata.obs['gRNA_count'] = [sgRNA_number_dict.get(x) for x in list(adata.obs_names)]
    adata.obs_names = [x+"_"+name for x in list(adata.obs_names)]
    adata.write(name+'.h5ad')
    return adata
    


# In[ ]:


gene_individual_files = []
for title in titles:
    original_path = os.path.join(base_dir, title, 'outs', 'filtered_feature_bc_matrix.h5')
    new_path = os.path.join(CRISPR_dir, f"{title}_filtered_feature_bc_matrix.h5")  # New file name based on title
    shutil.copy(original_path, new_path)
    gene_individual_files.append({'title': title, 'path': new_path})
h5ad_list = []   
for gene_file, guide_file in zip(gene_individual_files, individual_files):
    h5ad_list.append(process_data(gene_file['path'], guide_file['path'], ref_dict, gene_file['title']))


adata_merged = h5ad_list[0].concatenate(h5ad_list[1:], batch_key='sample', index_unique=None)
adata_merged.write("merged.h5ad")

