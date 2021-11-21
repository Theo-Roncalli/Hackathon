"""
Test stream pipeline to select the most variable
genes and run simulations departing from a raw
count matrix, instead of one that has been 
normalised and log-transformed.

In order to run this script you should have installed 
STREAM (https://github.com/pinellolab/STREAM)

Taken from STREAM's README.md :

$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge

$ conda create -n env_stream python stream=1.0 jupyter
$ conda activate env_stream

To install rna_seq_norm use :
pip install git+https://github.com/bnediction/RNA-seq-norm


Credits :
    " BNediction ;  pinellolab/STREAM "

Author:
    Gustavo Magaña López
"""

import numpy as np
import pandas as pd

import stream as st
import anndata as ann
import sys
import copy

# our development util :
from tagi_rnaseq.parsing import *
from tagi_rnaseq.normalization import *


def main(in_file, out_file, loess, percentile):
    """ intended for command line usage """
    # The whole raw dataset
    x = parse_counts_file(in_file)
    rename_counts_columns_in_place(x)

    # separate counts and metadata from the orignal dataframe
    counts = get_count_data(x)
    metadata = get_metadata(x)

    # normalize the counts
    cpm = normalize_cpm(counts)
    log_cpm = log_transform(cpm)
    y = log_cpm

    # use the STREAM pipeline to select the most variable genes :
    adata = ann.AnnData(y.T)
    adata.var_names_make_unique()
    st.set_workdir(adata, workdir="workdir")
    st.select_variable_genes(
        adata,
        loess_frac=loess,
        percentile=percentile,
        save_fig=True,
        fig_name=f"{in_file.split('/')[-1].split('.')[0]}_{loess}_{percentile}_std_vs_means.png",
    )

    # extract the genes index
    var_genes = adata.uns["var_genes"]

    # extract the raw counts
    raw_var = x.loc[var_genes, :]
    raw_var.to_csv(out_file)


if __name__ == "__main__":
    n_args = len(sys.argv) - 1
    if n_args == 4:
        main(sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]))
    else:
        print("")
        print("Select most variable genes.")
        print("usage: ")
        print(f"$ python {__file__} <<raw_counts.csv>> <<out_file.csv>> loess_frac percentile")
        print("example loess and percentile : ")
        print(f"$ python {__file__} <<raw_counts.csv>> <<out_file.csv>> 0.1 0.95")
        raise SystemExit("Incorrect usage.")
