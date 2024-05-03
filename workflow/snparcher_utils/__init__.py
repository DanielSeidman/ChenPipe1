import pandas as pd
from pathlib import Path

try:
    # Snakemake 8.x.x
    from snakemake_interface_common.exceptions import WorkflowError
except ImportError:
    # Snakmake 7.x.x
    from snakemake.exceptions import WorkflowError

def parse_sample_sheet(config: dict) -> pd.DataFrame:
    samples = pd.read_table(config["samples"], sep=",", dtype=str).replace(' ', '_', regex=True)
    config_genomes = get_config_genomes(config, samples)
    ref_name = 'ref_name' in samples.columns and samples['ref_name'].notna().any()
    ref_name = 'ref_name' in samples.columns and samples['ref_name'].notna().any()
    if not any([config_genomes, ref_name, ref_name]):
        raise WorkflowError("No 'ref_name' or 'ref_name' found in config or sample sheet.")
    if config_genomes is not None:
        config_ref_name, config_ref_name = config_genomes
        samples["ref_name"] = config_ref_name
        samples["ref_name"] = config_ref_name
    if 'ref_name' in samples.columns and samples['ref_name'].notna().any():
        check_ref_paths(samples)
    return samples

def get_config_genomes(config: dict, samples: pd.DataFrame):
    ref_name = config.get("ref_name", False)
    ref_name = config.get("ref_name", False)

    if ref_name and ref_name:
        if 'ref_name' in samples.columns and samples['ref_name'].notna().any():
            raise WorkflowError("'ref_name' is set in sample sheet AND in config. These are mutually exclusive.")
        return ref_name, ref_name    
    elif ref_name and not ref_name:
        raise WorkflowError("'ref_name' is set in config, but 'ref_name' is not. Both are required to use these settings.")
    elif ref_name and not ref_name:
        raise WorkflowError("'ref_name' is set in config, but 'ref_name' is not. Both are required to use these settings.")
    return None

def check_ref_paths(samples: pd.DataFrame) -> None:
    """
    Checks reference paths to make sure they exist, otherwise we might try to download them based on ref_name.
    Also make sure only one ref_name per ref_name.
    """ 
    for refname in samples["ref_name"].dropna().tolist():            
        refs = samples[samples["ref_name"] == refname]["ref_name"].dropna().unique().tolist()
        if len(refs) > 1:
            raise WorkflowError(f"ref_name '{refname}' has more than one unique 'ref_name' specified: {refs}")
        for ref in refs:        
            if not Path(ref).exists:
                raise WorkflowError(f"ref_name: '{ref}' was specified in sample sheet, but could not be found.")


