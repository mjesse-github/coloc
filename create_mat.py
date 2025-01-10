import os
import pandas as pd
import gzip
from tqdm import tqdm

def process_gwas_files_in_directory(base_directory, app_directory, summary_file_path):
    for root, dirs, files in os.walk(base_directory):
        for directory in dirs:
            dir_path = os.path.join(root, directory)
            for file in os.listdir(dir_path):
                if file.endswith('.lbf_variable.txt.gz') or file.endswith('.coloc5_combined.tsv'):
                    zip_flag = file.endswith('.lbf_variable.txt.gz')
                    process_gwas_file(os.path.join(dir_path, file), directory, app_directory, zip_flag, summary_file_path)

def process_gwas_file(file_path, directory_name, app_directory, zip_flag, summary_file_path):
    signals_dir = os.path.join(app_directory, "mat")
    os.makedirs(signals_dir, exist_ok=True)
    
    try:
        if zip_flag:
            with gzip.open(file_path, 'rt') as f:
                df = pd.read_csv(f, sep='\t', low_memory=False, on_bad_lines='skip')
        else:
            df = pd.read_csv(file_path, sep='\t', low_memory=False, on_bad_lines='skip')
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return

    trait_data = df.groupby(df['molecular_trait_id'])
    
    for trait_id, group in tqdm(trait_data, desc=f"Processing molecular traits for {file_path}"):
        for i in range(1, 11):
            lbf_column = f'lbf_variable{i}'
            if lbf_column in group.columns:
                process_signal(group, directory_name, trait_id, app_directory, lbf_column, i, summary_file_path)

def process_signal(group, directory_name, trait_id, app_directory, lbf_column, lbf_index, summary_file_path):
    df_filtered = group[['molecular_trait_id', 'region', 'variant', 'chromosome', 'position', lbf_column]].copy()
    df_filtered.rename(columns={lbf_column: 'lbf'}, inplace=True)
    
    if (df_filtered['lbf'] == 0).all():
        return
    
    signal_strength = df_filtered['lbf'].abs().max()
    if signal_strength < 1:
        return
    
    chromosome = df_filtered['chromosome'].iloc[0]
    location_min = df_filtered['position'].min()
    location_max = df_filtered['position'].max()

    signal = f"{directory_name}_{trait_id}_L{lbf_index}"
    output_file_name = f"{signal}.pickle"
    output_file_path = os.path.join(app_directory, output_file_name)

    df_filtered['lbf'] = pd.to_numeric(df_filtered['lbf'], errors='coerce')

    mat1_df = pd.DataFrame(df_filtered.set_index('variant')['lbf']).T

    mat1_df.to_pickle(output_file_path)

    summary_data = pd.DataFrame([{
        'signal': signal,
        'chromosome': chromosome,
        'location_min': location_min,
        'location_max': location_max,
        'signal_strength': signal_strength,
    }])
    
    header_needed = not os.path.exists(summary_file_path)  
    summary_data.to_csv(summary_file_path, sep='\t', mode='a', header=header_needed, index=False)

app_directory = "mat"
summary_file_path = "metadata.tsv"

base_directories = [
    "files"
]

for base_directory in base_directories:
    process_gwas_files_in_directory(base_directory, app_directory, summary_file_path)
