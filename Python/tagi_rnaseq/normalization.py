"""
docstring
"""

import numpy as np
import pandas as pd

def normalize_cpm(counts_df):
    """ Counts per million reads """
    return counts_df / counts_df.sum() * 1e6

def normalize_rpkm(counts_df, gene_lengths_series):
    """ """
    pass

def log_transform(data, base = None, constant = None):
    """
    log transform expression data

    Parameters
    ----------
    data: pd.DataFrame
        Annotated data matrix.
    base: int, optional (default: 2)
        The base used to calculate logarithm
    constant: int, optional (default: 1)
        The number to be added to all entries, to
        prevent indetermination of the transformation
        when a read is equal to zero.

    Returns
    -------
    a pandas.DataFrame, containing log-transformed counts
    """
    base = base or 2
    constant = constant or 1
    return np.log(data + constant) / np.log(base)

