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


# Setup 

You can run the setup by executing the management command `setup`:
```python
python manage.py setup
```

This will check if VEP, Typesense and MAGMA are ready for being used. If not, it will make sure that the necessary directories are created and the required files are downloaded.

# Preprocessing Steps

The preprocessing steps are necessary to prepare the data for the DyHealthNetLight backend.

## Normalize GWAS Summary Statistics

You can run this step by executing the management command `normalize`:
```python
python manage.py normalize
```

This creates the GWAS_stats_norm directory, which contains the normalized GWAS summary statistics.

## Generate Input Data: Manhattan, QQ, MAGMA files

You can run this step by executing the management command `input`:
```python
python manage.py input
```

This creates the GWAS_manhattan and GWAS_qq directories, which contain the generated Manhattan and QQ JSON files.

## VEP Annotation

To generate a VCF file containing all variants and annotate the variants using VEP, we plan to have the management command `annotate`. 

Once we implement our docker compose stack, this command will be changed.

You can run this step by executing the management command `annotate`:

```python
python manage.py annotate
```

## Typesense

You can run the typesense management command to fill the Typesense database with the necessary data:

```python
python manage.py typesense
```

## MAGMA Execution

### Download MAGMA

You can download the MAGMA software from the official website: https://ctg.cncr.nl/software/MAGMA

Make a directory called magma, download the zip and unzip it:

```bash
mkdir magma
cd magma
wget -O magma.zip https://vu.data.surfsara.nl/index.php/s/zkKbNeNOZAhFXZB/download
unzip magma.zip
```

Specify the location of the MAGMA executable in the .env file!

### Download Reference Data

You can download the reference data from the official website: https://ctg.cncr.nl/software/MAGMA#reference_data

Make a directory called magma_reference_data, download the zip and unzip it:

```bash
mkdir g1000_eur
cd g1000_eur
unzip g1000_eur.zip 
```

Specify the location of the LD reference data in the .env file!


# Run Django backend

After all preprocessing steps have been completed and the typesense docker container is up, you can run the Django backend.

```bash
python manage.py runserver 0.0.0.0:5136
```

Attention: the port needs to be the same as specified in the .env and remember to forward the ports.
