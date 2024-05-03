import glob
import re
import sys
import os
from pathlib import Path

def get_ref(wildcards):
    if "ref_name" in samples.columns:
        _refs = (
            samples.loc[(samples["ref_name"] == wildcards.ref_name)]["ref_name"]
            .dropna()
            .unique()
            .tolist()
        )
        if _refs:
            return _refs
    # if not user-specified ref_name, force MissingInputError in copy_ref with dummyfile, which allows download_ref to run b/c of ruleorder.
    logger.warning(f"snpArcher: ref_name specified in sample sheet header, but no path provided for ref_name '{wildcards.ref_name}'\n" +
                   f"Will try to download '{wildcards.ref_name}' from NCBI. If this is a genome accession, you can ignore this warning.")
    return []

# Get utils. This is not great, but we can move to setup.py and install via pip later if want
utils_path = (Path(workflow.main_snakefile).parent.parent.parent).resolve()
if str(utils_path) not in sys.path:
    sys.path.append(str(utils_path))

import pandas as pd
import snparcher_utils

samples = snparcher_utils.parse_sample_sheet(config)


