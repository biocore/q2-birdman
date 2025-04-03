# ----------------------------------------------------------------------------
# Copyright (c) 2024, Lucas Patel.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import pandas as pd
from joblib import Parallel, delayed
from qiime2 import Metadata
import biom
import numpy as np

from .src.birdman_chunked import run_birdman_chunk
from .src._utils import validate_table_and_metadata, validate_formula
from .src._summarize import summarize_inferences

def _create_dir(output_dir):
  sub_dirs = ["slurm_out", "logs", "inferences", "results", "plots"]
  for sub_dir in sub_dirs:
      os.makedirs(os.path.join(output_dir, sub_dir), exist_ok=True)

def run(table: biom.Table, metadata: Metadata, formula: str, threads: int = 16, 
        longitudinal: bool = False, subject_column: str = None) -> Metadata:
    """Run BIRDMAn and return the inference results as ImmutableMetadata."""
   
    validate_table_and_metadata(table, metadata)
    validate_formula(formula, table, metadata)
    
    metadata_df = metadata.to_dataframe()
    extra_params = {}
    
    # Only process longitudinal parameters if longitudinal=True
    if longitudinal:
        if subject_column is None:
            raise ValueError("subject_column must be specified when using longitudinal=True")
        if subject_column not in metadata_df.columns:
            raise ValueError(f"Subject column '{subject_column}' not found in metadata")
            
        group_var_series = metadata_df[subject_column]
        samp_subj_map = group_var_series.astype("category").cat.codes + 1
        groups = np.sort(group_var_series.unique())
        
        extra_params.update({
            "S": len(groups),
            "subj_ids": samp_subj_map.values,
            "u_p": 1.0  # Default value for subject random effects prior
        })
    
    chunks = 20
    output_dir = os.path.join(os.getcwd(), "test_out") 
    _create_dir(output_dir)
    print(f"Output directory is {output_dir}")
    os.makedirs(os.path.join(output_dir, "logs"), exist_ok=True)

    def run_chunk(chunk_num):
        log_path = os.path.join(output_dir, "logs", f"chunk_{chunk_num}.log")
        run_birdman_chunk(
            table=table,
            metadata=metadata_df,
            formula=formula,
            inference_dir=output_dir,
            num_chunks=chunks,
            chunk_num=chunk_num,
            logfile=log_path,
            longitudinal=longitudinal,
            **extra_params  # Only contains longitudinal params if longitudinal=True
        )

    Parallel(n_jobs=threads)(
        delayed(run_chunk)(i) for i in range(1, chunks + 1)
    )

    summarized_results = summarize_inferences(output_dir)
    summarized_results.index.name = 'featureid'
    results_metadata = Metadata(summarized_results)

    print(f"Results are stored in: {output_dir}")

    return results_metadata
