from pkg_resources import resource_filename
import biom
from birdman import SingleFeatureModel
import numpy as np
import pandas as pd

MODEL_PATH = resource_filename("q2_birdman.src", "stan/negative_binomial_lme_single.stan")

class ModelSingleLME(SingleFeatureModel):
    def __init__(
        self,
        table: biom.Table,
        feature_id: str,
        metadata: pd.DataFrame,
        formula: str,
        subj_ids: np.ndarray,
        S: int,
        beta_prior: float = 2.0,
        inv_disp_sd: float = 0.5,
        u_p: float = 1.0,
        absolute: bool = False,
        vi_iter=1000,
        num_draws=100,
        num_iter: int = 500,
        num_warmup: int = 500,
        **kwargs
    ):
        kwargs.pop('metadata', None)
        kwargs.pop('formula', None)
        kwargs.pop('absolute', None)

        super().__init__(
            table=table,
            feature_id=feature_id,
            model_path=MODEL_PATH,
            num_iter=num_iter,
            num_warmup=num_warmup,
            **kwargs
        )

        self.create_regression(formula=formula, metadata=metadata)

        if absolute:
            depth = np.zeros(table.shape[1])
            A = 0.0
        else:
            depth = np.log(table.sum(axis="sample"))
            A = np.log(1 / table.shape[0])

        param_dict = {
            "depth": depth,
            "B_p": beta_prior,
            "inv_disp_sd": inv_disp_sd,
            "A": A,
            "S": S,
            "subj_ids": subj_ids,
            "u_p": u_p
        }
        self.add_parameters(param_dict)

        self.specify_model(
            params=["beta_var", "inv_disp", "subj_int"],
            dims={
                "beta_var": ["covariate"],
                "subj_int": ["subject"],
                "log_lhood": ["tbl_sample"],
                "y_predict": ["tbl_sample"]
            },
            coords={
                "covariate": self.colnames,
                "subject": np.arange(1, S + 1),
                "tbl_sample": self.sample_names,
            },
            include_observed_data=True,
            posterior_predictive="y_predict",
            log_likelihood="log_lhood"
        )
