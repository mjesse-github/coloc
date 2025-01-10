#Jooksutab CLPP, filtreerib välja kõik, millel on rohkem

import gc
import gzip
import os

import pandas as pd

chromosome_column = 'chromosome' 
position_column = 'position'     
cs_id_column = 'cs_id'          
variant_column = 'variant'        

def process_directory(folder_path):
    dataframes = []
    
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith('credible_sets.tsv.gz'):
                filepath = os.path.join(root, file)
                with gzip.open(filepath, 'rt') as f:
                    df = pd.read_csv(f, sep='\t', low_memory=False, on_bad_lines='skip')
                if not df.empty:
                    df['dataset'] = file.split('.')[0]
                    dataframes.append(df)
    

    if dataframes:
        df_merged = pd.concat(dataframes, ignore_index=True)
        selected = df_merged[[cs_id_column, variant_column, 'pip', 'dataset']].rename(columns={'pip': 'pip_value'})
        selected.drop_duplicates(inplace=True)

        selected_a = selected.rename(columns={'pip_value': 'pip_A'})
        selected_a.drop(['dataset', cs_id_column], axis=1, inplace=True)
        
        selected_b = selected.rename(columns={'pip_value': 'pip_B'})
        selected_b.drop(['dataset', cs_id_column], axis=1, inplace=True)
        
        selected_a['signal_a'] = selected['dataset'] + "_" + selected[cs_id_column]
        selected_b['signal_b'] = selected['dataset'] + "_" + selected[cs_id_column]

        results = []
        for signal, group in selected_a.groupby("signal_a"):
            if group["variant"].nunique() < 400: 
                results.append(signal)

        selected_a = selected_a[selected_a["signal_a"].isin(results)]
        selected_b = selected_b[selected_b["signal_b"].isin(results)]

        selected_a.to_csv("selected.tsv", sep='\t', mode='a', header=True, index=False)
 
        output_file = 'cs_coloc.tsv'
        if os.path.exists(output_file):
            os.remove(output_file) 

        coloc2 = pd.merge(selected_a, selected_b, on='variant')

        coloc2 = coloc2.query('signal_a != signal_b')

        coloc2['pip_AB'] = coloc2['pip_A'] * coloc2['pip_B']
        coloc2 = coloc2.groupby(['signal_a', 'signal_b']).agg(clpp=('pip_AB', 'sum'), n_variants=('pip_AB', 'size')).reset_index()
        coloc2 = coloc2[coloc2['clpp'] > 0.04]

        coloc2.to_csv(output_file, sep='\t', mode='a', header=True, index=False)

        gc.collect()

folder_path = 'files' 
process_directory(folder_path)

