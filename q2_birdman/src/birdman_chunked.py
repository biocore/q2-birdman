import os
from tempfile import TemporaryDirectory
import time
import arviz as az
import biom
from birdman import ModelIterator
import cmdstanpy
import numpy as np
import pandas as pd
import logging
from .logger import setup_loggers
from .model_single import ModelSingle
from .model_single_lme import ModelSingleLME

logger = logging.getLogger(__name__)

def run_birdman_chunk(
    table,
    metadata,
    formula,
    inference_dir,
    num_chunks,
    chunk_num,
    chains=4,
    num_iter=500,
    num_warmup=500,
    beta_prior=5,
    inv_disp_sd=5,
    longitudinal=False,
    **kwargs
):
    FIDS = table.ids(axis="observation")
    birdman_logger = setup_loggers()

    # Choose appropriate model class and path
    if longitudinal:
        model_class = ModelSingleLME
    else:
        model_class = ModelSingle

    model_config = {
        "metadata": metadata,
        "formula": formula
    }

    model_kwargs = {
        "beta_prior": beta_prior,
        "inv_disp_sd": inv_disp_sd,
        "chains": chains,
        "num_iter": num_iter,
        "num_warmup": num_warmup
    }
    
    # Add longitudinal-specific parameters if needed
    if longitudinal:
        model_kwargs.update(kwargs)

    model_iter = ModelIterator(
        table,
        model_class,  # Use the appropriate model class
        num_chunks=num_chunks,
        **model_kwargs,
        **model_config
    )

    # Log the total number of chunks available
    total_chunks = len(model_iter)
    birdman_logger.info(f"ModelIterator created with {total_chunks} chunks.")

    # Check if the chunk_num is within range
    if chunk_num < 1 or chunk_num > total_chunks:
        birdman_logger.error(f"Chunk number {chunk_num} is out of range. Valid range is 1 to {total_chunks}.")
        return

    chunk = model_iter[chunk_num - 1]
    birdman_logger.info(f"Processing chunk number: {chunk_num}")

    for feature_id, model in chunk:
        feature_indices = np.where(FIDS == feature_id)[0]
        if len(feature_indices) == 0:
            birdman_logger.warning(f"Feature ID {feature_id} not found in FIDS")
            continue
        feature_num = feature_indices[0]
        feature_num_str = str(feature_num).zfill(4)
        birdman_logger.info(f"Processing feature {feature_id}")

        tmpdir = f"{inference_dir}/tmp/F{feature_num_str}_{feature_id}"
        infdir = f"{inference_dir}/inferences/"
        outfile = f"{inference_dir}/inferences/F{feature_num_str}_{feature_id}.nc"
        
        os.makedirs(infdir, exist_ok=True)
        os.makedirs(tmpdir, exist_ok=True)

        try:
            with TemporaryDirectory(dir=tmpdir) as t:
                try:
                    birdman_logger.info(f"Compiling model for feature {feature_id}")
                    model.compile_model()
                    birdman_logger.info(f"Fitting model for feature {feature_id}")
                    model.fit_model()
                    model.fit_model(sampler_args={"output_dir": t})
                except Exception as e:
                    birdman_logger.error(f"Error processing feature {feature_id}: {e}")
                    continue

                inf = model.to_inference()
                birdman_logger.info(f"Inference results for feature {feature_id}")

                # report diagnostics
                loo = az.loo(inf, pointwise=True)
                rhat = az.rhat(inf)
                if (rhat > 1.05).to_array().any().item():
                    birdman_logger.warning(f"{feature_id} has Rhat values > 1.05")
                if any(map(np.isnan, loo.values[:3])):
                    birdman_logger.warning(f"{feature_id} has NaN elpd values")

                # save inference results
                inf.to_netcdf(outfile)
                birdman_logger.info(f"Saved to {outfile}")
        except Exception as e:
            birdman_logger.error(f"Error with temporary directory for feature {feature_id}: {e}")
            continue
            
