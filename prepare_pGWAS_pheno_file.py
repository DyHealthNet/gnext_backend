import os
import pandas as pd


filepath = "/storage03/larend/pgwas_data/subtable0_20250730.txt"
df = pd.read_csv(filepath, sep="\t")
df["category"] = df["New.aptamer.in.7k.version"].map(
    lambda x: "new in 7k" if x else "previously assessed"
)

nr_samples = 13445

basepath = "/data/pheweb/pgwas/"

# loop over rows and check if the filename exists
new_df = pd.DataFrame(columns=["phenocode", "description", "nr_samples", "category", "external_id", "filename"])

for file in os.listdir(basepath):
    if file.endswith(".gwaslab.tsv.gz"):
        phenocode = file.split(".gwaslab.tsv.gz")[0]
        # subset df to seqid
        subset_df = df[df["SeqID"] == phenocode]
        description = " / ".join(subset_df["Target_Name"].astype(str).unique())
        external_id = " / ".join(subset_df["UniProt_ID"].astype(str).unique())
        category = " / ".join(subset_df["category"].astype(str).unique())

        # Create a one-row DataFrame to append
        new_row = pd.DataFrame([{
            "phenocode": phenocode,
            "description": description,
            "nr_samples": nr_samples,
            "category": category,
            "external_id": external_id,
            "filename": file
        }])

        # Append to existing DataFrame
        new_df = pd.concat([new_df, new_row], ignore_index=True)

new_df.to_csv("/storage03/larend/pgwas_data/pgwas_full_phenotype_file.csv", index=False)

