import logging
import os
import shutil
import struct
import tempfile
import time
import re
import lmdb
import msgpack
import pysam
from decouple import config
from backend.utils.preprocessing.zorp.zorp import parsers, sniffers, readers, lookups

logger = logging.getLogger("backend")

def map_and_write_rsid(norm_filepath, lmdb_path):

    if os.path.exists(norm_filepath):
        reader = sniffers.guess_gwas_standard(norm_filepath)
        status = add_rsid_to_gwas_stats(reader, norm_filepath, genome_build='GRCh37',
                                                   debug_mode=False, lmdb_path=lmdb_path)
    else:
        logger.info(f"Normalized GWAS stats file {norm_filepath} not available. Abort rsid mapping...")
        # TODO: exit?

def add_rsid_to_gwas_stats(reader, output_path, genome_build='GRCh37', debug_mode=False, lmdb_path=None):
    """
    Initial content ingestion: map rsID and write the new normalized file in standardized format with ending _rsid.
    If file is successfully written the original can be deleted if delete_original is set to true
    """
    delete_original = False  # Set to False if you want to keep the original

    build = lmdb_path + "/data.mdb" if lmdb_path else genome_build
    rsid_finder = lookups.SnpToRsid(build, test=debug_mode)
    reader.add_lookup('rsid', lambda variant: rsid_finder(variant.chrom, variant.pos, variant.ref, variant.alt))

    columns = ['chrom', 'pos', 'rsid', 'ref', 'alt', 'neg_log_pvalue', 'pvalue', 'beta', 'stderr_beta',
               'alt_allele_freq']

    i = 1
    for variant in reader:
        i+=1
        if i > 1:
            break

    # Build new output filename with _rsid suffix (.gz is not needed as write(...,make_tabix=True) will add that autom.)
    new_output_path = output_path.replace(".gz", "_rsid")
    if not os.path.exists(new_output_path + ".gz"):
        try:
            reader.write(new_output_path, make_tabix=True, columns=columns)
        except Exception as e:
            os.remove(new_output_path)
            logger.error(f"reader.write() failed. Deleted original file: {output_path}. Error: {e}")
            raise

    # Optionally delete original
    if delete_original:
        if os.path.exists(output_path):
            os.remove(output_path)
        if os.path.exists(output_path + ".tbi"):
            os.remove(output_path + ".tbi")
        logger.debug(f"Deleted original file and index: {output_path}")

def setup_rsid_mapping_lmdb(vcf_path, dir_path, num_chroms=25, map_size=10 ** 9):
    lmdb_path = dir_path + "/lmdb_" + config('VITE_GENOME_BUILD')
    # Only build the LMDB if it doesn't exist or is missing required files
    if not os.path.isdir(lmdb_path) or not os.path.exists(os.path.join(lmdb_path, "data.mdb")):
        logger.info("LMDB not found, creating...")
        start_time = time.time()
        build_snp_map_lmdb_from_vcf(vcf_path, lmdb_path)
        end_time = time.time()
        logger.debug(f"Time taken to produce LMDB mapping lib: {end_time - start_time:.2f} seconds")
    else:
        logger.info("LMDB already exists, skipping creation.")
    return lmdb_path

def build_snp_map_lmdb_from_vcf(vcf_path, lmdb_path, num_chroms=25):
    # Time taken to produce Lmdb mapping lib: 529.44 seconds/ 479.51 seconds
    # map_size is the maximum amount of data LMDB can store in this database., 1024 ** 4=1TB is a common upper limit
    env = lmdb.open(lmdb_path, map_size=1024 ** 4, max_dbs=num_chroms)
    db_handles = {}

    vcf = pysam.VariantFile(vcf_path)
    csq_header = vcf.header.info["CSQ"].description.split("Format: ")[1].split("|")
    rsid_index = csq_header.index("Existing_variation")

    with env.begin(write=True) as txn:
        for rec in vcf.fetch():
            chrom = rec.chrom

            # Open DB handle per chrom if not already open
            if chrom not in db_handles:
                db_handles[chrom] = env.open_db(chrom.encode(), txn=txn)

            db = db_handles[chrom]
            pos = rec.pos
            ref = rec.ref

            # Dictionary to hold ref/alt → rsID for this position
            refalt_to_rsid = {}

            for alt in rec.alts:

                # Parse rsID from CSQ
                csq_list = rec.info.get('CSQ')
                rsid = None
                if csq_list:
                    for csq_entry in csq_list:
                        fields = csq_entry.split('|')
                        if len(fields) <= rsid_index:
                            continue  # avoid index errors

                        # Get the Existing_variation field
                        candidate_field = fields[rsid_index]

                        # Split the candidate field on '&' to find valid rsIDs
                        for candidate in candidate_field.split('&'):
                            if candidate.startswith('rs') and candidate[2:].isdigit():
                                rsid = candidate
                                break
                        if rsid:
                            break

                if rsid is None:
                    # Use '.' or skip if no rsID found
                    continue

                refalt = f"{ref}/{alt}"
                # Store only the integer part of rsID, e.g. rs12345 → 12345
                rsid_int = int(rsid.replace('rs', ''))
                refalt_to_rsid[refalt] = rsid_int

            # If we found any ref/alt mappings for this position, store in LMDB
            if refalt_to_rsid:
                key_bytes = struct.pack('I', pos)  # pack position as 4-byte int
                value_bytes = msgpack.packb(refalt_to_rsid, use_bin_type=True)
                txn.put(key_bytes, value_bytes, db=db)

    env.sync()
    env.close()
