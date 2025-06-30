file_path = "/Users/lisiarend/Desktop/genes_pval/GS_UA_max_all_cohort_profile_paper_window_large2.71 .genes.out"

import pandas as pd

magma_output = pd.read_csv(file_path, sep='\s+')

# Filter significant genes based on p-value threshold
sig_thr = 0.05/len(magma_output)
print(len(magma_output))
print(sig_thr)

print(magma_output.columns)
significant_genes = magma_output[magma_output['P'] < sig_thr]

# Sort by p-value
significant_genes = significant_genes.sort_values(by='P')

# Print number of significant genes
print(f"Number of significant genes: {len(significant_genes)}")

print(significant_genes[['GENE', 'P']].head(10))
