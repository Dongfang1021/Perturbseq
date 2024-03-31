import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt


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

# Ensure the 'Report' directory exists
report_dir = os.path.join(base_dir, 'CRISPR')
os.makedirs(report_dir, exist_ok=True)

# Copy files into the 'Report' directory and generate a new list with updated paths
perturb_files = []
for title in titles:
    original_path = os.path.join(base_dir, title, 'outs', 'crispr_analysis', 'protospacer_calls_summary.csv')
    new_path = os.path.join(report_dir, f"{title}_protospacer_calls_summary.csv")  # New file name based on title
    shutil.copy(original_path, new_path)
    perturb_files.append({'title': title, 'path': new_path})
def basic_figure_total(file, name):    
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
	plt.savefig(name+'_feature_call.png')
       

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
	plt.savefig(name+'top10_feature_call.png')

	

for each in perturb_files:
	basic_figure_total(each['path'], each['title'])
	basic_figure_top10(each['path'], each['title'])
		