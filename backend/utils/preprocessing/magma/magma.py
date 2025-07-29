import csv
import re


def get_bool(value):
    val = value.strip().lower()
    true_values = {"true", "1", "yes", "on"}
    false_values = {"false", "0", "no", "off"}

    if val in true_values:
        return True
    elif val in false_values:
        return False
    else:
        return None

def sanitize_filename(s: str) -> str:
    """Sanitize string to be safe for filenames by replacing special chars. Equalise models that are synonymous in order
    to prevent redundant magma runs"""
    if "multi" in s:
        return "multi"
    if "snp-wise=top,1" == s:
        return "snp-wise=top"
    return re.sub(r'[^a-zA-Z0-9._-]', '_', s)

# MAGMA Configurations

# Allowed fixed model types
FIXED_ALLOWED_MODEL_TYPES = [
    "snp-wise=mean", "snp-wise=top", "multi", "multi=snp-wise", "snp-wise=multi"
]
# Allowed mapping strategies
ALLOWED_MAPPING_STRATEGIES = {"positional", "eQTL", "chromatin"} # or HiC

def is_valid_snpwise_top(key: str) -> bool:
    """Validate dynamic snp-wise=top,<value> where value is int > 0 or float between 0 and 1."""
    if not key.startswith("snp-wise=top,"):
        return False
    value = key.split("snp-wise=top,")[1]
    if value.isdigit() and int(value) > 0:
        return True
    try:
        float_val = float(value)
        return 0.0 < float_val <= 1.0
    except ValueError:
        return False

def is_allowed_model_type(key: str) -> bool:
    """Check if the model_type is allowed either fixed or dynamic snp-wise=top."""
    if key in FIXED_ALLOWED_MODEL_TYPES:
        return True
    return is_valid_snpwise_top(key)

def safe_non_negative_int(value, field_name, row):
   # if value in ('NA', 'Na', 'N/A', 'n/a', 'NULL', 'null', 'nan', 'NaN', 'missing'):
   #     return 0
    try:
        val = int(value)
        if val < 0:
            raise ValueError
        return val
    except (ValueError, TypeError):
        raise ValueError(
            f"Invalid value for '{field_name}' (must be a non-negative integer): {value} in row {row}"
        )

def read_magma_config(config_path="magma_config.csv", return_max_window=False):
    mconfigs = []
    seen = set()

    with open(config_path, newline="") as f:
        reader = csv.DictReader(f)
        window_up_max = 0
        window_down_max = 0
        for row in reader:
            model = row["model_type"].lower().strip()
            strategy = row["mapping_strategy"].lower().strip()

            # Parse window sizes, default 0 if missing or empty
            if strategy == "positional":
                window_up = safe_non_negative_int(row.get("window_up"), "window_up", row)
                window_down = safe_non_negative_int(row.get("window_down"), "window_down", row)
                if return_max_window:
                    window_up_max = window_up if window_up > window_up_max else window_up_max
                    window_down_max = window_down if window_down > window_down_max else window_down_max

            else:
                # Non-positional strategies don’t use windowing — default to 0 safely
                window_up = 0
                window_down = 0

            # Validate model_type
            if not is_allowed_model_type(model):
                raise ValueError(f"Unknown model_type '{model}' in config row: {row}")

            # Validate mapping_strategy
            if strategy not in ALLOWED_MAPPING_STRATEGIES:
                raise ValueError(f"Unknown mapping_strategy '{strategy}' in config row: {row}")

            # Skip identical rows (negligible)
            key = (model, strategy, window_up, window_down)
            if key in seen:
                continue
            seen.add(key)

            mconfigs.append({
                "model_type": model,
                "mapping_strategy": strategy,
                "window_up": window_up,
                "window_down": window_down
            })
    if return_max_window:
        return window_up_max, window_down_max

    return mconfigs

from decouple import Config, RepositoryEnv

# Example usage:
if __name__ == "__main__":
    config = Config(RepositoryEnv('/nfs/scratch/DyHealthNetLight/pheweb_backend_new/pheweb_backend/.env'))
    mconfigs = read_magma_config(config("MAGMA_CONFIG_FILE"))
    for cfg in mconfigs:
        print(cfg)