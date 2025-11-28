"""
LMDB Gene Query utilities for efficient variant-to-gene mapping.
"""

import logging
import struct
import msgpack
import lmdb

logger = logging.getLogger('backend')


class LMDBGeneQuery:
    """Context manager for efficient batch LMDB queries"""
    def __init__(self, lmdb_path):
        self.lmdb_path = lmdb_path
        self.env = None
    
    def __enter__(self):
        self.env = lmdb.open(self.lmdb_path, readonly=True, subdir=False, lock=False, max_dbs=128)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.env:
            self.env.close()
    
    def get_genes_for_variant(self, chrom, pos):
        """Query genes for a single variant. Environment stays open, but DB handles are not cached."""
        try:
            # Ensure chrom is a string and normalize
            chrom = str(chrom).replace('chr', '')
            
            # Ensure pos is an integer
            pos = int(pos)
            
            # Create key from position (4-byte unsigned big-endian integer)
            key_bytes = struct.pack(">I", pos)
            
            # Open database and query in same transaction (don't cache handles)
            with self.env.begin(write=False) as txn:
                try:
                    db = self.env.open_db(chrom.encode(), txn=txn)
                    value_bytes = txn.get(key_bytes, db=db)
                except lmdb.NotFoundError:
                    logger.warning(f"Chromosome database '{chrom}' not found in LMDB")
                    return []
            
            # Unpack result
            if value_bytes:
                genes = msgpack.unpackb(value_bytes, raw=False)
                # Sort by distance (already sorted, but explicit for clarity)
                genes_sorted = sorted(genes, key=lambda x: x[2])
                # Return only symbols
                return {ensg_id:symbol for ensg_id, symbol, distance in genes_sorted}
            return {}
        except lmdb.Error as e:
            logger.error(f"LMDB error for {chrom}:{pos} - {e}")
            return {}
        except Exception as e:
            logger.error(f"Error getting genes for variant {chrom}:{pos} - {e}")
            return {}

