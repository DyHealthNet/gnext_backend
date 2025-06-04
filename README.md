# Backend of the PheWeb Extension Platform

# Conda Environment

To create the conda environment, run the following command:

```bash 
conda env create -f environment.yml
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
 docker run -d --name typesense-server -p 8108:8108 -v pheweb_typesense_data:/data typesense/typesense:29.0.rc15 --data-dir /data --api-key=xyz --enable-cors
```

Run the typesense_initialization.py script to fill the database with the initial data (currently: Pan_UKBB)
```bash
python typesense_initialization.py
```


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

Annotate file with VEP:

```bash
docker run -v /nfs/data/Pan_UKBB/data:/data ensemblorg/ensembl-vep vep --cache --offline --input_file --output_file -y GRCh37 
```