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
import logging

from .src.birdman_chunked import run_birdman_chunk
from .src._utils import validate_table_and_metadata, validate_formula
from .src._summarize import summarize_inferences

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def _create_dir(output_dir):
    sub_dirs = ["slurm_out", "logs", "inferences", "results", "plots"]
    for sub_dir in sub_dirs:
        os.makedirs(os.path.join(output_dir, sub_dir), exist_ok=True)

def run(table: biom.Table, metadata: Metadata, formula: str, threads: int = 16, 
        longitudinal: bool = False, subject_column: str = None,
        beta_prior: float = 5.0, inv_disp_sd: float = 5.0, u_p: float = 1.0) -> Metadata:
    """Run BIRDMAn and return the inference results as ImmutableMetadata."""
   
    validate_table_and_metadata(table, metadata)
    validate_formula(formula, table, metadata)
    
    metadata_df = metadata.to_dataframe()
    extra_params = {
        "beta_prior": beta_prior,
        "inv_disp_sd": inv_disp_sd
    }
    
    # Ensure number of threads doesn't exceed number of features
    num_features = table.shape[0]
    if threads > num_features:
        logger.warning(f"Number of threads ({threads}) exceeds number of features ({num_features}). "
              f"Reducing threads to {num_features}.")
        threads = num_features
    
    if longitudinal:
        if subject_column is None or subject_column not in metadata_df.columns:
            raise ValueError(f"Subject column '{subject_column}' must be specified and exist in metadata when using longitudinal=True")
            
        group_var_series = metadata_df[subject_column]
        samp_subj_map = group_var_series.astype("category").cat.codes + 1
        groups = np.sort(group_var_series.unique())
        
        extra_params.update({
            "S": len(groups),
            "subj_ids": samp_subj_map.values,
            "u_p": u_p  # subject random effects prior
        })
    
    chunks = min(20, num_features)
    logger.info(f"Processing {num_features} features in {chunks} chunks using {threads} threads")
    
    with tempfile.TemporaryDirectory() as output_dir:
        _create_dir(output_dir)
        logger.info(f"Using temporary directory: {output_dir}")

        def run_chunk(chunk_num):
            logger.info(f"Starting chunk {chunk_num}")
            run_birdman_chunk(
                table=table,
                metadata=metadata_df,
                formula=formula,
                inference_dir=output_dir,
                num_chunks=chunks,
                chunk_num=chunk_num,
                longitudinal=longitudinal,
                **extra_params
            )
            logger.info(f"Completed chunk {chunk_num}")

        # Use loky backend for better process management
        results = Parallel(n_jobs=threads, backend='loky', verbose=10)(
            delayed(run_chunk)(i) for i in range(1, chunks + 1)
        )
        logger.info("All chunks completed")

        summarized_results = summarize_inferences(output_dir)
        summarized_results.index.name = 'featureid'
        results_metadata = Metadata(summarized_results)

        return results_metadata
