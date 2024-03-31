from jinja2 import Environment, FileSystemLoader
from bs4 import BeautifulSoup
import os
import shutil

# Define the base directory where the subfolders are located
base_dir = '/run/media/Dongfang/DataUSB/Perturb_seq/10XGenomixs_10/'

# Define the list of HTML files based on titles, assuming each has its own subfolder named after the title
titles = [
    'NoTNFa_CS1_U6d1_A3', 'NoTNFa_CS1_10X_B3', 'TNFa_U6d1_B2_75',
    'TNFa_U6d1_B2_25', 'TNFa_CS1_U6d1_B4', 'NoTNFa_U6d1_A1',
    'TNFa_CS1_10X_B4', 'NoTNFa_CS1_10X_A3', 'TNFa_CS1_10X_A4',
    'NoTNFa_CS1_U6d1_B3', 'TNFa_U6d1_A2', 'TNFa_CS1_U6d1_A4',
    'NoTNFa_U6d1_B1_25', 'NoTNFa_U6d1_B1_75'
]

# Ensure the 'Report' directory exists
report_dir = os.path.join(base_dir, 'Report')
os.makedirs(report_dir, exist_ok=True)

# Copy files into the 'Report' directory and generate a new list with updated paths
html_files = []
for title in titles:
    original_path = os.path.join(base_dir, title, 'outs', 'web_summary.html')
    new_path = os.path.join(report_dir, f"{title}.html")  # New file name based on title
    shutil.copy(original_path, new_path)
    html_files.append({'title': title, 'path': new_path})
    
combined_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Combined HTML Report</title>
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css">
</head>
<body>
    <div class="container mt-3">
        <ul class="nav nav-tabs" id="myTab" role="tablist">"""

for index, file in enumerate(html_files):
    tab_id = f"tab{index}"
    combined_html += f"""
            <li class="nav-item">
                <a class="nav-link{' active' if index == 0 else ''}" id="{tab_id}-tab" data-toggle="tab" href="#{tab_id}" role="tab">{file['title']}</a>
            </li>"""

combined_html += """
        </ul>
        <div class="tab-content" id="myTabContent">"""

for index, file in enumerate(html_files):
    tab_id = f"tab{index}"
    combined_html += f"""
            <div class="tab-pane fade{' show active' if index == 0 else ''}" id="{tab_id}" role="tabpanel" aria-labelledby="{tab_id}-tab">
                <iframe src="{os.path.basename(file['path'])}" style="width:100%;height:100vh;border:none;"></iframe>
            </div>"""

combined_html += """
        </div>
    </div>
    <script src="https://code.jquery.com/jquery-3.5.1.slim.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/@popperjs/core@2.5.2/dist/umd/popper.min.js"></script>
    <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/js/bootstrap.min.js"></script>
</body>
</html>"""

# Write the combined HTML content to a new file
with open(os.path.join(report_dir, 'combined_report.html'), 'w') as file:
    file.write(combined_html)

print("Combined HTML report generated in the 'Report' directory.")