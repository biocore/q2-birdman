import pandas as pd
import numpy as np
from qiime2 import Metadata
from q2_types.metadata import ImmutableMetadata
from q2_types.feature_data import FeatureData, Taxonomy
from .plugin_setup import plugin

@plugin.register_transformer
def _1(data: ImmutableMetadata) -> Metadata:
    """
    Transform an ImmutableMetadata artifact to a Metadata object.
    
    Parameters
    ----------
    data : ImmutableMetadata
        The ImmutableMetadata artifact to transform
        
    Returns
    -------
    Metadata
        A Metadata object containing the data from the ImmutableMetadata artifact
    """
    return data.view(Metadata)

@plugin.register_transformer
def _2(data: FeatureData[Taxonomy]) -> pd.DataFrame:
    """
    Transform a FeatureData[Taxonomy] artifact to a pandas DataFrame.
    
    Parameters
    ----------
    data : FeatureData[Taxonomy]
        The taxonomy artifact to transform
        
    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the taxonomy data
    """
    return data.view(pd.DataFrame)