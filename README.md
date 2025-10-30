# Backend of the DyHealthNetLight Platform

# Conda Environment

Create a conda environment using the requirements.txt file

```bash
conda create -n dyhealthnetlight_backend python=3.10
conda activate dyhealthnetlight_backend
pip install -r requirements.txt
```

Activate your conda environment

```bash
conda activate dyhealthnetlight_backend
```


# Preprocessing Pipeline

Before running the platform, you need to preprocess your GWAS summary statistics file using our Nextflow pipeline, available here: https://github.com/DyHealthNet/dyhealthnetlight_nf_pipeline


# Setup

## Configuration

You need to fill the .env file with your study-specific parameters. For checking whether everything in the .env file is correct, you can run the check_env management command.

```python
python manage.py check_env
```

## Typesense

Before being able to deploy the platform, you can run the typesense management command to download the Typesense docker image, run the container, and fill the Typesense database with the necessary data:

```python
python manage.py typesense
```


# Run Django backend

After all preprocessing steps have been completed and the typesense docker container is up, you can run the Django backend.

```bash
python manage.py runserver 0.0.0.0:5136
```

Attention: the port needs to be the same as specified in the .env and remember to forward the ports.
