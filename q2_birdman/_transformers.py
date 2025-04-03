import pandas as pd
import numpy as np
from qiime2 import Metadata
from q2_types.metadata import ImmutableMetadata
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