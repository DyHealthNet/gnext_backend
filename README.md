# Backend of the PheWeb Extension Platform

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



# Typesense Initialization

First, pull the docker image:

```bash
docker pull typesense/typesense:29.0.rc15  
```

Then, create a volume:

```bash
docker volume create pheweb_typesense_data
```

Finally, run the container:

```bash
 docker run -d --name typesense-server -p 8108:8108 -v pheweb_typesense_data:/data typesense/typesense:29.0.rc15 --data-dir /data --backend-key=xyz --enable-cors
```

Run the typesense_initialization.py script to fill the database with the initial data (currently: Pan_UKBB)
```bash
python typesense_initialization.py
```

# Generation of VCF file with all variants

For Pan_UKBB, the Bash script `make_variant_vcf.sh` (in `/nfs/data/Pan_UKBB/`) can be used to generate a VCF file with all variants.

# Annotation of Variants with VEP

Make a vep_data directory:
```bash
mkdir vep_data
```

Pull the docker image:
```bash
docker pull ensemblorg/ensembl-vep
```

Download the VEP cache:

Ensure that Docker has access to your data directory.

```bash
docker run -t -i -v /nfs/data/Pan_UKBB/data:/data ensemblorg/ensembl-vep INSTALL.pl -a cf -s homo_sapiens -y GRCh37
```

Cache stored in --fasta /opt/vep/.vep/homo_sapiens/113_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz

For the annotation, the Bash script `vep_run.sh` (in `/nfs/data/Pan_UKBB/`) can be used. This creates the `annotated_variants.vcf.bgz` file.


# Generation of Manhattan and QQ JSON Files

For this, the api/utils/locuszoom/reader.py file needs to be run. This creates new directories with files per GWAS summary statistic.

# Run Django backend

After all preprocessing steps have been completed and the typesense docker container is up, you can run the Django backend.

```bash
python manage.py runserver 0.0.0.0:5136
```

Attention: the port needs to be the same than specified in the .env and remember to forward the ports.
