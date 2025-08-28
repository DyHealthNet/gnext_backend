# build_chr_alltraits_bgz_rowblocks_parallel.py
import os
import re
import argparse
import logging
import numpy as np
import pandas as pd
import pysam
from math import ceil
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import shared_memory
from decouple import config
from django.conf import settings

logger = logging.getLogger("backend")

# 0-based columns in per-trait Tabix files:
# CHR(0) POS(1) ID(2) REF(3) ALT(4) neg_log_pvalue(5) BETA(6) SE(7) AF(8)
COLS = {
    "neg_log_pvalue": 5,
    "beta": 7,
    "stderr_beta": 8,
    "alt_allele_freq": 9,
}
METRICS = ["neg_log_pvalue", "beta", "stderr_beta", "alt_allele_freq"]

def load_traits_and_paths():
    """
        Load trait names and corresponding file paths from the phenotype CSV.

        Reads the file specified by the PHENO_FILE environment variable,
        extracts the 'phenocode' column as trait names, and constructs the absolute paths
        to the normalized GWAS result files in GWAS_NORM_DIR.

        Returns:
            tuple:
                traits (list[str]): List of trait codes as strings.
                paths (list[str]): List of absolute file paths to per-trait Tabix-indexed GWAS files.
    """
    pheno = pd.read_csv(config("PHENO_FILE"))
    traits = pheno["phenocode"].astype(str).tolist()
    paths = [
        os.path.join(
            settings.GWAS_NORM_DIR,
            re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(i)) + ".gz"
        )
        for i in pheno["filename"].astype(str).tolist()
    ]
    return traits, paths

def build_chr_variant_index(chrom: str):
    """
        Build a variant index for a specific chromosome from the master annotated VCF.

        Opens the GWAS_ANNO_VCF_FILE in GWAS_VEP_DIR using pysam.VariantFile, iterates
        through all records for the specified chromosome, and stores variant identifiers
        and coordinates.

        Args:
            chrom (str): Chromosome name as in the VCF (e.g., "1" or "chr1").

        Returns:
            tuple:
                rows (list[tuple]): List of (chrom, pos, ref, alt, vid_bytes) for each variant.
                row_index (dict): Mapping from vid_bytes to row index in 'rows'.
    """
    vcf_path = os.path.join(settings.GWAS_VEP_DIR, settings.GWAS_ANNO_VCF_FILE)
    rows = []
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf.fetch(str(chrom)):
            c = str(rec.contig)
            p = int(rec.pos)
            r = rec.ref
            a = rec.alts[0]
            vid = f"{c}:{p}:{r}:{a}".encode()
            rows.append((c, p, r, a, vid))
    row_index = {vid: i for i, (_, _, _, _, vid) in enumerate(rows)}
    return rows, row_index

def _parse_float(s: str) -> float:
    """
        Parse a string into a float, returning NaN on failure.

        Args:
            s (str): Input string to convert.

        Returns:
            float: Parsed float value, or numpy.nan if parsing fails.
    """
    try:
        return float(s)
    except Exception:
        return np.nan

def open_writers(chr_out, trait_names):
    """
        Open BGZF-compressed TSV writers for all metrics for a chromosome.

        Creates the output directory for the chromosome, writes a common header
        line containing 'chrom', 'pos', 'ref', 'alt', 'vid', and all trait names
        to each metric file.

        Args:
            chr_out (str): Path to chromosome-specific output directory.
            trait_names (list[str]): List of trait names.

        Returns:
            tuple:
                paths (dict): Mapping from metric name to output file path.
                handles (dict): Mapping from metric name to open BGZF file handle.
    """
    os.makedirs(chr_out, exist_ok=True)
    header = "#" + "\t".join(["chrom", "pos", "ref", "alt", "vid"] + trait_names) + "\n"
    paths = {
        "neg_log_pvalue": os.path.join(chr_out, "neg_log_pvalue.tsv.bgz"),
        "beta": os.path.join(chr_out, "beta.tsv.bgz"),
        "stderr_beta": os.path.join(chr_out, "stderr_beta.tsv.bgz"),
        "alt_allele_freq": os.path.join(chr_out, "alt_allele_freq.tsv.bgz"),
    }
    handles = {m: pysam.BGZFile(paths[m], "w") for m in paths}
    for h in handles.values():
        h.write(header.encode())
    return paths, handles

def close_and_index(paths, handles):
    """
        Close BGZF writers and create Tabix indexes for each metric file.

        Args:
            paths (dict): Mapping from metric name to output file path.
            handles (dict): Mapping from metric name to open BGZF file handle.
        """
    for m, h in handles.items():
        h.close()
        pysam.tabix_index(paths[m], seq_col=0, start_col=1, end_col=1, meta_char="#", force=True)

def _worker_fill_traits(chrom, trait_idxs, trait_paths, row_index, i0, i1,
                        shm_names, shape, dtype, cols):
    """
        Worker function to fill metric values for a subset of traits in a row block.

        Attaches to shared memory arrays for each metric, opens each trait's Tabix
        file, fetches all variants for the given chromosome, and fills metric arrays
        for rows in the range [i0, i1).

        Args:
            chrom (str): Chromosome name.
            trait_idxs (list[int]): List of trait indices to process in this worker.
            trait_paths (list[str]): Absolute paths to all trait files.
            row_index (dict): Mapping from vid_bytes to row index.
            i0 (int): Start row index for the block (inclusive).
            i1 (int): End row index for the block (exclusive).
            shm_names (dict): Shared memory names for each metric array.
            shape (tuple): Shape of each metric array (rows, traits).
            dtype (numpy.dtype): Data type of metric arrays.
            cols (tuple[int]): Column indices for (pval, beta, stderr, af).
    """
    import os as _os
    import numpy as _np
    import pysam as _pysam
    col_p, col_b, col_se, col_af = cols

    p_shm = shared_memory.SharedMemory(name=shm_names['p'])
    b_shm = shared_memory.SharedMemory(name=shm_names['b'])
    se_shm = shared_memory.SharedMemory(name=shm_names['se'])
    af_shm = shared_memory.SharedMemory(name=shm_names['af'])

    p_arr = _np.ndarray(shape, dtype=dtype, buffer=p_shm.buf)
    b_arr = _np.ndarray(shape, dtype=dtype, buffer=b_shm.buf)
    se_arr = _np.ndarray(shape, dtype=dtype, buffer=se_shm.buf)
    af_arr = _np.ndarray(shape, dtype=dtype, buffer=af_shm.buf)

    for j in trait_idxs:
        fn = trait_paths[j]
        if not (_os.path.exists(fn) and _os.path.exists(fn + ".tbi")):
            continue
        try:
            with _pysam.TabixFile(fn) as tbx:
                try:
                    it = tbx.fetch(str(chrom))
                except ValueError:
                    continue
                for line in it:
                    f = line.rstrip("\n").split("\t")
                    if len(f) <= max(col_p, col_b, col_se, col_af, 4):
                        continue
                    vid = f"{f[0]}:{f[1]}:{f[3]}:{f[4]}".encode()
                    gi = row_index.get(vid)
                    if gi is None or gi < i0 or gi >= i1:
                        continue
                    li = gi - i0
                    # fast parse
                    try:
                        p_arr[li, j] = float(f[col_p])
                    except Exception:
                        p_arr[li, j] = _np.nan
                    try:
                        b_arr[li, j] = float(f[col_b])
                    except Exception:
                        b_arr[li, j] = _np.nan
                    try:
                        se_arr[li, j] = float(f[col_se])
                    except Exception:
                        se_arr[li, j] = _np.nan
                    try:
                        af_arr[li, j] = float(f[col_af])
                    except Exception:
                        af_arr[li, j] = _np.nan
        except Exception:
            continue
        logger.info(f"[Worker] chrom {chrom} trait with index {j} proceesed in shard")

    p_shm.close(); b_shm.close(); se_shm.close(); af_shm.close()

def fill_block_parallel(chrom, i0, i1, trait_paths, row_index,
                        p_block, b_block, se_block, af_block,
                        max_workers):
    """
        Fill a row block's metric arrays in parallel across multiple trait shards.

        Creates shared memory segments for the four metric arrays, splits traits
        into contiguous shards, launches _worker_fill_traits() for each shard, and
        copies data back from shared memory into the local arrays.

        Args:
            chrom (str): Chromosome name.
            i0 (int): Start row index for the block (inclusive).
            i1 (int): End row index for the block (exclusive).
            trait_paths (list[str]): Absolute paths to all trait files.
            row_index (dict): Mapping from vid_bytes to row index.
            p_block (numpy.ndarray): Array for p-values.
            b_block (numpy.ndarray): Array for beta values.
            se_block (numpy.ndarray): Array for standard errors.
            af_block (numpy.ndarray): Array for allele frequencies.
            max_workers (int): Number of parallel workers to use.
    """
    ntraits = p_block.shape[1]
    block_rows = i1 - i0
    dtype = np.float32

    def _mk_shm(arr):
        shm = shared_memory.SharedMemory(create=True, size=arr.nbytes)
        np.ndarray(arr.shape, dtype=arr.dtype, buffer=shm.buf)[:] = arr
        return shm

    p_shm = _mk_shm(p_block)
    b_shm = _mk_shm(b_block)
    se_shm = _mk_shm(se_block)
    af_shm = _mk_shm(af_block)
    shm_names = {'p': p_shm.name, 'b': b_shm.name, 'se': se_shm.name, 'af': af_shm.name}

    # contiguous shards across traits
    shards = []
    if max_workers <= 1:
        shards = [list(range(ntraits))]
    else:
        step = ceil(ntraits / max_workers)
        for s in range(max_workers):
            j0 = s * step
            j1 = min((s + 1) * step, ntraits)
            if j0 < j1:
                shards.append(list(range(j0, j1)))

    cols = (COLS["neg_log_pvalue"], COLS["beta"], COLS["stderr_beta"], COLS["alt_allele_freq"])

    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futs = [
            ex.submit(
                _worker_fill_traits, chrom, shard, trait_paths, row_index, i0, i1,
                shm_names, (block_rows, ntraits), dtype, cols
            )
            for shard in shards
        ]
        for f in as_completed(futs):
            f.result()

    # copy shared back to local arrays
    p_arr = np.ndarray((block_rows, ntraits), dtype=dtype, buffer=p_shm.buf)
    b_arr = np.ndarray((block_rows, ntraits), dtype=dtype, buffer=b_shm.buf)
    se_arr = np.ndarray((block_rows, ntraits), dtype=dtype, buffer=se_shm.buf)
    af_arr = np.ndarray((block_rows, ntraits), dtype=dtype, buffer=af_shm.buf)
    p_block[:] = p_arr
    b_block[:] = b_arr
    se_block[:] = se_arr
    af_block[:] = af_arr

    p_shm.close(); p_shm.unlink()
    b_shm.close(); b_shm.unlink()
    se_shm.close(); se_shm.unlink()
    af_shm.close(); af_shm.unlink()

def build_chromosome(chrom: str, out_dir: str, row_block_size: int, trait_workers: int):
    """
        Build BGZF+Tabix metric files for all traits for a given chromosome.

        Loads trait names and paths, builds the chromosome variant index, opens
        BGZF writers for each metric, processes variants in row blocks, fills
        block arrays in parallel across traits, and writes rows to the output
        files. Closes and indexes the metric files when done.

        Args:
            chrom (str): Chromosome name (as in VCF).
            out_dir (str): Output base directory.
            row_block_size (int): Number of variants to hold in memory per block.
            trait_workers (int): Number of parallel workers to use for trait processing.
    """
    trait_names, trait_paths = load_traits_and_paths()
    ntraits = len(trait_names)

    rows, row_index = build_chr_variant_index(chrom)
    nvar = len(rows)
    if nvar == 0:
        logger.info(f"{chrom}: no variants, skip")
        return

    chr_out = os.path.join(out_dir, f"chr{chrom}")
    paths, writers = open_writers(chr_out, trait_names)

    logger.info(f"chr{chrom}: variants={nvar}, traits={ntraits}, row_block_size={row_block_size}, trait_workers={trait_workers}")

    missing = [fn for fn in trait_paths if not (os.path.exists(fn) and os.path.exists(fn + ".tbi"))]
    if missing:
        logger.warning(f"{len(missing)} trait files missing bgz or .tbi, they will be skipped")

    nblocks = ceil(nvar / row_block_size)
    for b in range(nblocks):
        i0 = b * row_block_size
        i1 = min((b + 1) * row_block_size, nvar)
        block_rows = i1 - i0

        p_block = np.full((block_rows, ntraits), np.nan, dtype=np.float32)
        b_block = np.full((block_rows, ntraits), np.nan, dtype=np.float32)
        se_block = np.full((block_rows, ntraits), np.nan, dtype=np.float32)
        af_block = np.full((block_rows, ntraits), np.nan, dtype=np.float32)

        fill_block_parallel(
            chrom, i0, i1, trait_paths, row_index,
            p_block, b_block, se_block, af_block,
            max_workers=max(1, trait_workers)
        )

        for li in range(block_rows):
            c, p, r, a, vid = rows[i0 + li]
            head = [c, str(p), r, a, vid.decode()]
            writers["neg_log_pvalue"].write(("\t".join(head + [("." if not np.isfinite(x) else f"{x:.6g}") for x in p_block[li, :]]) + "\n").encode())
            writers["beta"].write(("\t".join(head + [("." if not np.isfinite(x) else f"{x:.6g}") for x in b_block[li, :]]) + "\n").encode())
            writers["stderr_beta"].write(("\t".join(head + [("." if not np.isfinite(x) else f"{x:.6g}") for x in se_block[li, :]]) + "\n").encode())
            writers["alt_allele_freq"].write(("\t".join(head + [("." if not np.isfinite(x) else f"{x:.6g}") for x in af_block[li, :]]) + "\n").encode())

        logger.info(f"chr{chrom}: wrote rows [{i0}:{i1})")

    close_and_index(paths, writers)
    logger.info(f"chr{chrom}: done")
