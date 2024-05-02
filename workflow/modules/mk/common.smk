import glob
import re
import sys
import os
from pathlib import Path

# Get utils. This is not great, but we can move to setup.py and install via pip later if want
utils_path = (Path(workflow.main_snakefile).parent.parent.parent).resolve()
if str(utils_path) not in sys.path:
    sys.path.append(str(utils_path))

import pandas as pd
import snparcher_utils

samples = snparcher_utils.parse_sample_sheet(config)


