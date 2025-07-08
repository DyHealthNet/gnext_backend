# Backend of the DyHealthNetLight Platform

# Conda Environment

Create a conda environment using the requirements.txt file

```bash
conda create -n pheweb_backend python=3.10
conda activate pheweb_backend
pip install -r requirements.txt
```

Activate your conda environment

```bash
conda activate pheweb_backend
```

# Preprocessing Steps

The preprocessing steps are necessary to prepare the data for the DyHealthNetLight backend.

# Normalize GWAS Summary Statistics

You can run this step by executing the management command `normalize`:
```python
python manage.py normalize
```

This creates the GWAS_stats_norm directory, which contains the normalized GWAS summary statistics.

# Generate Manhattan and QQ Plots

You can run this step by executing the management command `locuszoom`:
```python
python manage.py locuszoom
```

This creates the GWAS_manhattan and GWAS_qq directories, which contain the generated Manhattan and QQ JSON files.

# VEP Annotation

To generate a VCF file containing all variants and annotate the variants using VEP, we plan to have the management command `annotate`. 

Since this command requires docker, we currently run this step manually. Once we implement our docker compose stack, this command will be available.

Currently, you need to run the following steps manually:

## Create a VCF file with all variants

```bash
mkdir GWAS_vep_directory
backend/utils/preprocessing/bash/generate_full_variants_vcf.sh GWAS_norm_dir GWAS_vcf_file
```

## Setup VEP

```bash
backend/utils/preprocessing/bash/setup_vep.sh GWAS_vep_directory genome_build
```

## Annotate variants with VEP

```bash
backend/utils/preprocessing/bash/run_vep.sh GWAS_vcf_file GWAS_annotated_vcf_file GWAS_vcf_dir genome_build
```


# Typesense Initialization

To initialize the typesense engine for the autocomplete functionality in the frontend, we plan to have the management command `typesense`. 

Since this command requires docker, we currently run this step manually. Once we implement our docker compose stack, this command will be available.

Currently, you need to run the following steps manually:

First, pull the docker image:

```bash
docker pull typesense/typesense:29.0.rc15  
```

Then, create a volume:

```bash
docker volume create GWAS_typesense
```

Finally, run the container:

```bash
 docker run -d --name typesense-server -p 8108:8108 -v GWAS_typesense:/data typesense/typesense:29.0.rc15 --data-dir /data --backend --api-key=xyz --enable-cors
```

Now you can run the typesense management command to initialize the typesense engine:

```python
python manage.py typesense
```

# Run Django backend

After all preprocessing steps have been completed and the typesense docker container is up, you can run the Django backend.

```bash
python manage.py runserver 0.0.0.0:5136
```

Attention: the port needs to be the same as specified in the .env and remember to forward the ports.

# MAGMA Execution

## Download MAGMA

Magma gets downloaded automatically. If, however, the download fails or is outdated, it is also possible to download it yourself and specify the path to the executable in the .env.

You can download the MAGMA software from the official website: https://ctg.cncr.nl/software/MAGMA

Make a directory called magma, download the zip and unzip it:

```bash
mkdir magma
cd magma
wget -O magma.zip https://vu.data.surfsara.nl/index.php/s/zkKbNeNOZAhFXZB/download
unzip magma.zip
```

Specify the location of the MAGMA executable in the .env file!

## Download Reference Data

You can download the reference data from the official website: https://ctg.cncr.nl/software/MAGMA#reference_data

Make a directory called magma_reference_data, download the zip and unzip it:

```bash
mkdir g1000_eur
cd g1000_eur
unzip g1000_eur.zip 
```

Specify the location of the LD reference data in the .env file!

## Run MAGMA

TODO