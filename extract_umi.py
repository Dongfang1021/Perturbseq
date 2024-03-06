"""
extract mean umi from single cell data

Dongfang Hu

dfhu@cellecta.com

20240304


"""

import scanpy as sc
import argparse
import os
from datetime import datetime
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--sc_input', help='h5ad files contain sc data')
parser.add_argument('--output', help='output file data path and file name')
parser.add_argument('--gene_list', help='gene list used to extract mean umi for each gene')
parser.add_argument('--subset', help='subset of sc data used to extract mean umi for each gene in provided gene list')
argv = vars(parser.parse_args())

if argv['sc_input'] == None:
    	raise Exception ('You should provide sc h5ad!')
else:
	sc_input = argv['sc_input'].strip()

if argv['output'] == None:
    	raise Exception ('You should provide output file name')
else:
	output=argv['output'].strip()

if argv['gene_list'] == None:
    	raise Exception ('You should provide a gene list!')
else:
	gene_list=argv['gene_list'].strip()




#TODO add args input output
#TODO save files as a table
#TODO extract sub dataset default is for all the cells







# Load your .h5ad file
adata = sc.read_h5ad(sc_input)
# load gene list
gene_list_df = pd.read_csv(gene_list, sep='\t')
# Define your list of genes of interest
genes_of_interest = list(gene_list_df['transcript_symbol'])  # Replace these with your actual genes

# Check if the genes are in the dataset
missing_genes = [gene for gene in genes_of_interest if gene not in adata.var_names]
if missing_genes:
    print("Warning: The following genes are missing from the dataset:", missing_genes)

# Calculate mean UMI counts for the genes of interest
mean_umi_counts = {}
for gene in genes_of_interest:
    if gene in adata.var_names:
        gene_index = list(adata.var_names).index(gene)
        mean_umi_counts[gene] = adata.X[:, gene_index].mean()

# Convert dictionary to DataFrame
df = pd.DataFrame(mean_umi_counts, index=[0])
df.to_csv(output, sep='\t', index=False)
