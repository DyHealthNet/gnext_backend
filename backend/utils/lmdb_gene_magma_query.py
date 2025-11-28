"""
LMDB Gene MAGMA Query utilities for efficient gene MAGMA p-value lookups.
"""

import logging
import lmdb
import pysam
import gzip

logger = logging.getLogger('backend')


class LMDBGeneMAGMAQuery:
    """Context manager for efficient gene MAGMA p-value lookups by ENSG_ID using line numbers"""
    
    def __init__(self, lmdb_path, bgz_path):
        """
        Initialize the query manager.
        
        Args:
            lmdb_path: Path to LMDB index file mapping ENSG_ID to line number
            bgz_path: Path to BGZ file containing gene MAGMA p-values
        """
        self.lmdb_path = lmdb_path
        self.bgz_path = bgz_path
        self.env = None
        self.trait_names = None
    
    def __enter__(self):
        # Open LMDB environment
        self.env = lmdb.open(self.lmdb_path, readonly=True, subdir=False, lock=False)
        
        # Read header to get trait names (only once)
        with pysam.BGZFile(self.bgz_path, "r") as f:
            header_line = f.readline().decode('utf-8').strip()
            if header_line.startswith('#'):
                header_line = header_line.lstrip('#')
            columns = header_line.split('\t')
            self.trait_names = columns[2:]  # Skip GENE and SYMBOL
        
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.env:
            self.env.close()
    
    def get_gene_magma_pvalues(self, ensg_id):
        """
        Query MAGMA p-values for a single gene using LMDB-indexed virtual offset.
        Uses pysam's seek() for O(1) random access to BGZ file.
        
        Args:
            ensg_id: Ensembl gene ID (e.g., "ENSG00000123456")
            
        Returns:
            dict with keys:
                - 'ensg_id': Gene ID
                - 'symbol': Gene symbol
                - 'trait_pvalues': list of dicts with 'trait_id' and 'pvalue'
            Returns None if gene not found
        """
        try:
            # Query LMDB for virtual offset
            with self.env.begin(write=False) as txn:
                offset_bytes = txn.get(ensg_id.encode('utf-8'))
            
            if offset_bytes is None:
                logger.warning(f"Gene '{ensg_id}' not found in MAGMA index")
                return None
            
            # Convert bytes to integer (8-byte little-endian)
            virtual_offset = int.from_bytes(offset_bytes, byteorder='little')
            
            # Use pysam's seek() to jump directly to the gene's line (O(1) random access!)
            with pysam.BGZFile(self.bgz_path, "r") as f:
                f.seek(virtual_offset)
                line = f.readline()
                
                if not line:
                    logger.error(f"Could not read at virtual offset {virtual_offset} for gene {ensg_id}")
                    return None
                
                line = line.decode('utf-8').strip()
                
                # Parse line
                columns = line.split('\t')
                gene_id = columns[0]
                symbol = columns[1]
                pvalues_raw = columns[2:]
                
                # Verify gene ID matches
                if gene_id != ensg_id:
                    logger.error(f"Gene ID mismatch: expected {ensg_id}, got {gene_id}")
                    return None
                
                # Build trait p-value list
                trait_pvalues = []
                for trait_id, pval_str in zip(self.trait_names, pvalues_raw):
                    if pval_str == '.':
                        pvalue = None
                    else:
                        try:
                            pvalue = float(pval_str)
                        except ValueError:
                            pvalue = None
                    
                    trait_pvalues.append({
                        'trait_id': trait_id,
                        'pvalue': pvalue
                    })
                
                return {
                    'ensg_id': gene_id,
                    'symbol': symbol,
                    'trait_pvalues': trait_pvalues
                }
            
        except lmdb.Error as e:
            logger.error(f"LMDB error for gene {ensg_id}: {e}")
            return None
        except Exception as e:
            logger.error(f"Error getting MAGMA p-values for gene {ensg_id}: {e}")
            return None
    
    def get_trait_names(self):
        """
        Get list of all trait IDs in the MAGMA file.
        
        Returns:
            list of trait ID strings
        """
        return self.trait_names.copy() if self.trait_names else []
