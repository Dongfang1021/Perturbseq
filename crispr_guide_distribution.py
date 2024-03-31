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
report_dir = os.path.join(base_dir, 'CRISPR')
os.makedirs(report_dir, exist_ok=True)

# Copy files into the 'Report' directory and generate a new list with updated paths
perturb_files = []
for title in titles:
    original_path = os.path.join(base_dir, title, 'outs', 'crispr_analysis', 'protospacer_calls_per_cell.csv')
    new_path = os.path.join(report_dir, f"{title}_protospacer_calls_per_cell.csv")  # New file name based on title
    shutil.copy(original_path, new_path)
    perturb_files.append({'title': title, 'path': new_path})

reference = pd.read_csv('RDL103_feature_ref.csv')
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
	plt.savefig(name+'_sgRNA_number_distribution.png')
       
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
	plt.savefig(name+'_umi_distribution.png', bbox_inches='tight')
	plt.cla()

	plt.figure(figsize=(24, 6))
	sns.boxplot(x='feature_call', y='num_umis', data=detected_single_guide_umi, palette='bright')
	# Save the figure with bbox_inches='tight' to include all labels
	plt.xticks(rotation=90, fontsize=8)
	plt.yscale('log')
	plt.title(name+'UMI Distribution.png')
	plt.savefig(name+'_umi_distribution_boxplot.png', bbox_inches='tight')
	plt.cla()


	plt.figure(figsize=(24, 6))
	sns.lineplot(x='feature_call', y='mean umi', data=detected_single_guide_umi_mean,  palette='bright', size=3)
	plt.yscale('log')
	plt.gca().invert_yaxis()

	# Save the figure with bbox_inches='tight' to include all labels
	plt.xticks(rotation=90, fontsize=8)
	plt.title(name+'Mean UMI Distribution.png')
	plt.savefig(name+'_mean_umi_distribution.png', bbox_inches='tight')
	return len(guide_set - detected_single_guide_set)

missing_guide_dict = {}	

for each in perturb_files:
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
plt.savefig('Missing_guide_number.png')

		