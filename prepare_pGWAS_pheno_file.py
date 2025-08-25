import os
import pandas as pd


filepath = "/storage03/larend/pgwas_data/subtable0_20250730.txt"
df = pd.read_csv(filepath, sep="\t")

#print(df.head())

#print(df.columns)


basepath = "/data/pheweb/pgwas/"

# loop over rows and check if the filename exists
new_df = pd.DataFrame(columns=["phenocode", "description", "category", "external_id", "filename"])

for file in os.listdir(basepath):
    if file.endswith(".gwaslab.tsv.gz"):
        phenocode = file.split(".gwaslab.tsv.gz")[0]
        # subset df to seqid
        subset_df = df[df["SeqID"] == phenocode]
        category = " / ".join(subset_df["class"].astype(str).unique())
        description = " / ".join(subset_df["HARMONIZED_GENE_NAME"].astype(str).unique())
        external_id = " / ".join(subset_df["UniProt_ID"].astype(str).unique())

        if category == "nan":
            category = "Not given"

        if description == "nan":
            description = "Not given"

        # Create a one-row DataFrame to append
        new_row = pd.DataFrame([{
            "phenocode": phenocode,
            "description": description,
            "category": category,
            "external_id": external_id,
            "filename": file
        }])

        # Append to existing DataFrame
        new_df = pd.concat([new_df, new_row], ignore_index=True)

new_df.to_csv("/storage03/larend/pgwas_data/pgwas_full_phenotype_file.csv", index=False)

