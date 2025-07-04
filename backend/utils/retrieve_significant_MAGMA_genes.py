#file_path = "/Users/lisiarend/Desktop/genes_pval/GS_UA_max_all_cohort_profile_paper_window.genes.out"
file_path = "/Users/lisiarend/Downloads/magma.genes.outfiles/UA_magma.genes.out"

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

significant_genes.to_csv("/Users/lisiarend/Desktop/UA_significant_genes.csv", index=False)


# FDR control

#from statsmodels.stats.multitest import multipletests

#rejected, pvals_corrected, _, _ = multipletests(magma_output['P'], alpha=0.05, method='fdr_bh')
#magma_output['FDR'] = pvals_corrected
#significant_genes = magma_output[magma_output['FDR'] < 0.05]

# Sort by p-value
#significant_genes = significant_genes.sort_values(by='FDR')

# Print number of significant genes
#print(f"Number of significant genes: {len(significant_genes)}")

#print(significant_genes[['GENE', 'P', 'FDR']].head(10))
