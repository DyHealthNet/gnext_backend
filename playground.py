import pandas as pd
from pathlib import Path
import os

# Path to directory and manifest
dir_path = Path("/storage03/larend/pgwas_data/new_pgwas_out/GWAS_manhattan")
manifest = pd.read_csv("/storage03/larend/pgwas_data/pgwas_full_phenotype_file.csv")
all_phenos = manifest["phenocode"].tolist()

# Extract existing phenocodes from filenames
existing = [p.name.split(".gwaslab_manhattan")[0] for p in dir_path.glob("*.gwaslab_manhattan.json")]

# Compare
missing = sorted(set(all_phenos) - set(existing))

for m in missing:
    print(m)
print(f"Missing ({len(missing)}):")

missing = sorted(set(existing) - set(all_phenos))
print(f"Missing ({len(missing)}):")
