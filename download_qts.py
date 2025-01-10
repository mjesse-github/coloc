import os
import pandas as pd
import urllib.request
import time

metadata_path = 'tabix_ftp_paths.tsv'  
metadata = pd.read_csv(metadata_path, sep='\t')  

metadata = metadata[(metadata["quant_method"].isin(["ge", "microarray"])) & (metadata["study_label"]=="GTEx")]

download_dir = 'files'
os.makedirs(download_dir, exist_ok=True)

def download_and_unzip_file(url, gz_path):
    while True: 
        try:
            urllib.request.urlretrieve(url, gz_path)
            print(f"Downloaded {gz_path}")
            break 
        except Exception as e:
            print(f"Error downloading {url}: {e}. Retrying in 5 seconds...")
            time.sleep(5)  

for index, row in metadata.iterrows():
    study_id = row['study_id']
    dataset_id = row['dataset_id']

    print(f"Starting {dataset_id}")

    os.makedirs(f"{download_dir}/{dataset_id}", exist_ok=True)

    gz_path = os.path.join(download_dir,dataset_id, f"{dataset_id}.credible_sets.tsv.gz")
    tsv_path = os.path.join(download_dir,dataset_id, f"{dataset_id}.credible_sets.tsv")
    
    download_and_unzip_file(row["ftp_cs_path"], gz_path)

    gz_path = os.path.join(download_dir,dataset_id, f"{dataset_id}.lbf_variable.txt.gz")
    tsv_path = os.path.join(download_dir,dataset_id, f"{dataset_id}.lbf_variable.txt.gz")
    download_and_unzip_file(row["ftp_lbf_path"], gz_path)

