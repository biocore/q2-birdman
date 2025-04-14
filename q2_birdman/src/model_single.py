from pkg_resources import resource_filename

import biom
# IF USING A DIFFERENT STAN MODEL, CHANGE HERE
from birdman import SingleFeatureModel
import numpy as np
import pandas as pd

# PROVIDE FILEPATH TO STAN MODEL
# ENSURE THAT YOU HAVE COMPILED THIS MODEL BY 
# ACTIVATING BIRDMAN ENVIRONMENT
# OPENING PYTHON/IPYTHON  
# IMPORTING CMDSTANPY, THEN RUNNING 
# cmdstanpy.CmdStanModel(stan_file="path/to/model.stan")
MODEL_PATH = resource_filename("q2_birdman.src", "stan/negative_binomial_single.stan")

# NAME CLASS SOMETHING RELEVANT TO YOUR MODEL
class ModelSingle(SingleFeatureModel):
    def __init__(
        self,
        table: biom.Table,
        feature_id: str,
        metadata: pd.DataFrame,
        formula: str,
        # OPTIONAL: CHANGE PARAMETERS
        beta_prior: float = 2.0,
        inv_disp_sd: float = 0.5,
        vi_iter=1000,
        num_draws=100,
        num_iter: int = 500,
        num_warmup: int = 500,
        **kwargs
    ):

        kwargs.pop('metadata', None)  # remove metadata if it's in kwargs, not sure
        kwargs.pop('formula', None)   # remove formula if it's in kwargs, not sure

        super().__init__(
            table=table,
            feature_id=feature_id,
            model_path=MODEL_PATH,
            num_iter=num_iter,
            num_warmup=num_warmup,
            **kwargs
        )

        self.create_regression(formula=formula, metadata=metadata)

        param_dict = {
            "depth": np.log(table.sum(axis="sample")),
            "B_p": beta_prior,
            "inv_disp_sd": inv_disp_sd,
	          "A": np.log(1 / table.shape[0])
        }
        self.add_parameters(param_dict)

        self.specify_model(
            params=["beta_var", "inv_disp"],
            dims={
                "beta_var": ["covariate"],
                "log_lhood": ["tbl_sample"],
                "y_predict": ["tbl_sample"]
            },
            coords={
                "covariate": self.colnames,
                "tbl_sample": self.sample_names,
            },
            include_observed_data=True,
            posterior_predictive="y_predict",
            log_likelihood="log_lhood"
        )
