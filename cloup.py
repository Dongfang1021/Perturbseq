import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import numpy as np
import seaborn as sns

# Define the base directory where the subfolders are located
base_dir = '/run/media/Dongfang/DataUSB/Perturb_seq/10XGenomixs_10/'

# Define the list of HTML files based on titles, assuming each has its own subfolder named after the title
# titles = [
#     'NoTNFa_CS1_U6d1_A3', 'NoTNFa_CS1_10X_B3', 'TNFa_U6d1_B2_75',
#     'TNFa_U6d1_B2_25', 'TNFa_CS1_U6d1_B4', 'NoTNFa_U6d1_A1',
#     'TNFa_CS1_10X_B4', 'NoTNFa_CS1_10X_A3', 'TNFa_CS1_10X_A4',
#     'NoTNFa_CS1_U6d1_B3', 'TNFa_U6d1_A2', 'TNFa_CS1_U6d1_A4',
#     'NoTNFa_U6d1_B1_25'
# ]
titles = [
    'NoTNFa_CS1_U6d1_A3', 'NoTNFa_CS1_10X_B3', 'TNFa_U6d1_B2_75',
    'TNFa_U6d1_B2_25', 'TNFa_CS1_U6d1_B4', 'NoTNFa_U6d1_A1',
    'TNFa_CS1_10X_B4', 'NoTNFa_CS1_10X_A3', 'TNFa_CS1_10X_A4',
    'NoTNFa_CS1_U6d1_B3', 'TNFa_U6d1_A2', 'TNFa_CS1_U6d1_A4',
    'NoTNFa_U6d1_B1_25','NoTNFa_U6d1_B1_75'
]
titles = sorted(titles)
# Ensure the 'Report' directory exists
report_dir = os.path.join(base_dir, 'Cloupe')
os.makedirs(report_dir, exist_ok=True)

# Copy files into the 'Report' directory and generate a new list with updated paths
perturb_files = []
for title in titles:
    original_path = os.path.join(base_dir, title, 'outs',  'cloupe.cloupe')
    new_path = os.path.join(report_dir, f"{title}_cloupe.cloupe")  # New file name based on title
    shutil.copy(original_path, new_path)
    perturb_files.append({'title': title, 'path': new_path})