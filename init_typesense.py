#!/usr/bin/env python3
"""
Standalone script to initialize Typesense collections and fill with data.
This runs AFTER Typesense container is started but BEFORE backend serves requests.
"""
import sys
import os
import time
import traceback
import logging
import json
import gzip
import requests
import pandas as pd
import typesense
from decouple import config

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("init_typesense")


def wait_for_typesense(host: str, port: str, api_key: str, timeout: int = 600):
    """Wait for Typesense to become healthy (may take time when loading large collections from disk)."""
    url = f"http://{host}:{port}/health"
    headers = {"X-TYPESENSE-API-KEY": api_key}

    logger.info(f"Waiting for Typesense at {url} (timeout: {timeout}s)...")
    start = time.time()
    while time.time() - start < timeout:
        try:
            r = requests.get(url, headers=headers, timeout=5)
            if r.ok:
                logger.info("✅ Typesense is ready.")
                return True
        except requests.exceptions.RequestException:
            pass
        time.sleep(2)
    raise TimeoutError("Typesense did not become ready in time.")


def create_schema(client):
    """Create or verify autocomplete collection schema."""
    schema_autocomplete = {
        "name": "autocomplete",
        "fields": [
            {"name": "id", "type": "string"},
            {"name": "type", "type": "string", "facet": True},
            {"name": "label", "type": "string"},
            {"name": "description", "type": "string"},
            {"name": "external_ref", "type": "string"},
            {"name": "category", "type": "string"},
            {"name": "filename", "type": "string"},
            {"name": "nr_samples", "type": "int32"},
        ],
    }

    # Create collection if missing
    collections = client.collections.retrieve()
    if "autocomplete" not in [collection["name"] for collection in collections]:
        try:
            client.collections.create(schema_autocomplete)
            logger.info("Typesense autocomplete schema created!")
        except Exception as e:
            raise RuntimeError(f"Failed to create Typesense schema: {e}")
    else:
        logger.info("Typesense autocomplete schema already exists, skipping creation.")


def reset_collection_if_needed(client, collection_name, force_reset=False):
    """Reset collection if it has existing data (or skip if data exists and force_reset=False)."""
    try:
        results = client.collections[collection_name].documents.search({
            "q": "*",
            "query_by": "description",
            "per_page": 1
        })

        if results["found"] > 0:
            if not force_reset:
                logger.info(f"Collection {collection_name} has {results['found']} documents. Skipping re-import (set FORCE_TYPESENSE_RESET=true to force).")
                return False  # Signal to skip import
            
            logger.info(f"Collection {collection_name} has {results['found']} documents. Force resetting...")
            client.collections[collection_name].delete()
            schema_autocomplete = {
                "name": "autocomplete",
                "fields": [
                    {"name": "id", "type": "string"},
                    {"name": "type", "type": "string", "facet": True},
                    {"name": "label", "type": "string"},
                    {"name": "description", "type": "string"},
                    {"name": "external_ref", "type": "string"},
                    {"name": "category", "type": "string"},
                    {"name": "filename", "type": "string"},
                    {"name": "nr_samples", "type": "int32"},
                ],
            }
            client.collections.create(schema_autocomplete)
            logger.info(f"Collection {collection_name} reset successfully.")
        else:
            logger.info(f"No existing documents found in {collection_name}.")
        
        return True  # Signal to proceed with import
    except Exception as e:
        logger.error(f"Error resetting collection {collection_name}: {e}")
        return True  # Proceed with import on error


def import_phenotypes(client, pheno_file):
    """Import phenotypes to Typesense."""
    logger.info(f"Importing phenotypes from {pheno_file}")
    pheno_dt = pd.read_csv(pheno_file)
    
    for i, r in pheno_dt.iterrows():
        logger.info(f"Importing phenotype to typesense: {r['phenocode']}")
        if "external_id" not in r:
            r["external_id"] = ""
        if "nr_samples" not in r:
            r["nr_samples"] = 0

        doc = {
            "type": "trait",
            "id": str(r["phenocode"]),
            "label": str(r["phenocode"]),
            "description": r["description"],
            "external_ref": r["external_id"],
            "category": r["category"],
            "filename": r["filename"],
            "nr_samples": int(r["nr_samples"]),
        }
        client.collections["autocomplete"].documents.upsert(doc)
    
    logger.info(f"Imported {len(pheno_dt)} phenotypes.")


def import_genes(client, gene_file):
    """Import genes to Typesense."""
    logger.info(f"Reading gene file: {gene_file}")
    mapped_genes = pd.read_csv(gene_file, sep="\t")
    logger.info(f"Loaded {len(mapped_genes)} genes")
    
    documents = []
    for i, r in mapped_genes.iterrows():
        doc = {
            "type": "gene",
            "id": str(r["ensg_id"]),
            "label": str(r["symbol"]),
            "description": f"Chr {r['chr']}: {r['start']}-{r['end']} ({r['strand']})",
            "external_ref": str(r["symbol"]),
            "category": "",
            "filename": "",
            "nr_samples": 0,
        }
        documents.append(doc)
    
    logger.info(f"Prepared {len(documents)} gene documents for import")
    payload = "\n".join(json.dumps(doc) for doc in documents)
    result = client.collections["autocomplete"].documents.import_(payload, {"action": "upsert"})
    logger.info(f"Genes imported to typesense! Result: {result}")


def import_variants(client, vcf_file, batch_size):
    """Import variants to Typesense."""
    logger.info(f"Importing variants from {vcf_file}")
    
    # First, parse the annotation columns from the VCF header
    anno_columns = []
    with gzip.open(vcf_file, "rt") as vcf:
        for line in vcf:
            if line.startswith("##INFO") and "Format:" in line:
                anno_columns = line.strip().split("|")
                anno_columns[0] = anno_columns[0].split("Format: ")[1].split('"')[0]
                anno_columns[-1] = anno_columns[-1].strip('">')
                break

    # Now process the variants
    documents = []
    batch_nr = 1
    header = []

    with gzip.open(vcf_file, "rt") as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                header = [h.replace("#", "") for h in line.strip().split("\t")]
                continue
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            variant_dict = dict(zip(header, fields))

            # Parse INFO field into structured dictionaries
            info = variant_dict["INFO"].replace("CSQ=", "").split(",")
            info_dict_list = [dict(zip(anno_columns, i.split("|"))) for i in info]

            # Extract unique rs_ids from the parsed dictionaries
            rs_ids = list(
                set([d.get("Existing_variation", "") for d in info_dict_list if d.get("Existing_variation", "")]))
            rs_ids = ", ".join(rs_ids)

            doc = {
                "type": "variant",
                "id": f"{variant_dict['CHROM']}_{variant_dict['POS']}_{variant_dict['REF']}/{variant_dict['ALT']}",
                "label": f"{variant_dict['CHROM']}_{variant_dict['POS']}_{variant_dict['REF']}/{variant_dict['ALT']}",
                "description": f"{variant_dict['CHROM']}_{variant_dict['POS']}_{variant_dict['REF']}/{variant_dict['ALT']}",
                "external_ref": rs_ids,
                "category": "",
                "filename": "",
                "nr_samples": 0,
            }
            documents.append(doc)

            if len(documents) % batch_size == 0:
                payload = "\n".join(json.dumps(doc) for doc in documents)
                client.collections["autocomplete"].documents.import_(payload, {"action": "upsert"})
                logger.info(f"Batch Nr. {batch_nr} imported to typesense!")
                documents.clear()
                batch_nr += 1

        # Import remaining documents
        if documents:
            payload = "\n".join(json.dumps(doc) for doc in documents)
            client.collections["autocomplete"].documents.import_(payload, {"action": "upsert"})
            logger.info(f"Batch Nr. {batch_nr} imported to typesense!")


def main():
    try:
        # Get configuration from environment
        TYPESENSE_HOST = config("VITE_TYPESENSE_HOST")
        TYPESENSE_PORT = config("VITE_TYPESENSE_PORT")
        TYPESENSE_KEY = config("VITE_TYPESENSE_KEY")
        PHENO_FILE = config("PHENO_FILE")
        BATCH_SIZE = int(config("BATCH_SIZE"))
        
        # Get file paths from Django-like settings
        NF_DATA_DIR = config("NF_DATA_DIR")
        ANNO_VCF_FILE = os.path.join(NF_DATA_DIR, "annotate/full_variants.vcf.gz")
        GENE_FILE = os.path.join(NF_DATA_DIR, "lmdb_gene/mapped_genes.tsv")
        
        logger.info("=" * 70)
        logger.info("Typesense Initialization Script")
        logger.info("=" * 70)
        logger.info(f"Host: {TYPESENSE_HOST}:{TYPESENSE_PORT}")
        logger.info(f"Phenotype file: {PHENO_FILE}")
        logger.info(f"Gene file: {GENE_FILE}")
        logger.info(f"VCF file: {ANNO_VCF_FILE}")
        logger.info("=" * 70)
        
        # Wait for Typesense to be ready
        wait_for_typesense(TYPESENSE_HOST, TYPESENSE_PORT, TYPESENSE_KEY)
        
        # Initialize client
        client = typesense.Client({
            "nodes": [{
                "host": TYPESENSE_HOST,
                "port": TYPESENSE_PORT,
                "protocol": "http"
            }],
            "api_key": TYPESENSE_KEY,
            "connection_timeout_seconds": 200
        })
        
        # Create schema
        create_schema(client)
        
        # Check if we should reset/import data
        force_reset = config("FORCE_TYPESENSE_RESET", default="false").lower() == "true"
        should_import = reset_collection_if_needed(client, "autocomplete", force_reset=force_reset)
        
        if should_import:
            # Import data
            import_phenotypes(client, PHENO_FILE)
            import_genes(client, GENE_FILE)
            import_variants(client, ANNO_VCF_FILE, BATCH_SIZE)
        else:
            logger.info("Skipping data import - collection already has data.")
        
        logger.info("=" * 70)
        logger.info("✅ Typesense initialization completed successfully!")
        logger.info("=" * 70)
        
    except Exception as e:
        logger.error("=" * 70)
        logger.error("❌ Typesense initialization failed!")
        logger.error("=" * 70)
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
