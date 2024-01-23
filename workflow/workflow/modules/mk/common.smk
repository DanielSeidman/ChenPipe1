import glob
import re
import os
import sys
import pandas as pd
from snakemake.exceptions import WorkflowError

samples = pd.read_table(config["samples"], sep=",", dtype=str).replace(' ', '_', regex=True)


