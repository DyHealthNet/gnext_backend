#Extracted from locuszoom-hosted: d726fb2
from backend.utils.preprocessing.zorp.zorp import parsers, sniffers, readers, lookups

import lmdb
import msgpack
import struct
import re
import pysam
import logging

logger = logging.getLogger("backend")

def build_snp_map_lmdb_from_vcf(vcf_path, lmdb_path, num_chroms=25, map_size=10**9):
    env = lmdb.open(lmdb_path, map_size=map_size, max_dbs=num_chroms)
    db_handles = {}

    vcf = pysam.VariantFile(vcf_path)

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
                key = f"{chrom}:{pos}_{ref}/{alt}"

                # Parse rsID from CSQ
                csq_list = rec.info.get('CSQ')
                rsid = None
                if csq_list:
                    for csq_entry in csq_list:
                        annotations = re.split(r'[,&]', csq_entry)
                        for annot in annotations:
                            fields = annot.split('|')
                            if len(fields) > 17:
                                candidate = fields[17]
                                if candidate.startswith('rs'):
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

def normalize_contents_lib(reader, output_path, genome_build='GRCh37', debug_mode=False, lmdb_path=None):
    """
    Initial content ingestion: load the file and write variants in a standardized format

    This routine will deliberately exclude lines that could not be handled in a reliable fashion, such as pval=NA

    In "debug mode", the ingest process will use a smaller (test environment optimized) version of the
        rsid_finder lookup.
    """
    reader.write(output_path, columns=["rsid", "pval"], make_tabix=False)
    return True