# yourapp/management/commands/build_chr_bgz.py
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

# import your builder
# if you used the row-block parallel version:
from backend.utils.preprocessing.build_chr_metric_bgz import build_chromosome

class Command(BaseCommand):
    help = "Build per-chromosome BGZF files for all traits and 4 metrics"

    def add_arguments(self, parser):
        parser.add_argument("--chrom", required=True, help="Chromosome as in VCF, e.g. 1 or chr1")
        parser.add_argument("--out-dir", default=getattr(settings, "GWAS_CHR_BGZ_DIR", "./bgz_out"))
        parser.add_argument("--row-block-size", type=int, default=64000)
        parser.add_argument("--trait-workers", type=int, default=8)

    def handle(self, *args, **opts):
        chrom = opts["chrom"]
        out_dir = opts["out_dir"] or getattr(settings, "GWAS_CHR_BGZ_DIR", "./bgz_out")
        row_block_size = int(opts["row_block_size"])
        trait_workers = int(opts["trait_workers"])

        # Call the builder
        build_chromosome(
            chrom=chrom,
            out_dir=out_dir,
            row_block_size=row_block_size,
            trait_workers=trait_workers,
        )