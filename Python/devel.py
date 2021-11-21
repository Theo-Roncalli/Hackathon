
try:
    from profile_binr import ProfileBin
except:
    print("\n ImportError (catched) : unable to import profile_binr")

from tagi_rnaseq.parsing import *
from tagi_rnaseq.normalization import *

# read counts
x = parse_counts_file("../Data/counts.txt")
rename_counts_columns_in_place(x)

# separate counts and metadata from the orignal dataframe
counts = get_count_data(x)
metadata = get_metadata(x)

# normalize the counts
cpm = normalize_cpm(counts)
log_cpm = log_transform(cpm)

# most variable genes
vargenes = pd.read_csv("var_counts.csv", index_col=0)
