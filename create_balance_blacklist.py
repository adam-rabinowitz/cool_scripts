import argparse
import cooler
import numpy as np


def get_clr_counts(clr):
    '''Calculates sums and non-zero pixels for each row/col'''
    # Get coo and get shape
    coo = clr.matrix(sparse=True, balance=False)[:]
    nrow, ncol = coo.shape
    assert(nrow == ncol)
    # Check blacklist:
    if blacklist_pixels:
        assert(min(blacklist_pixels) >= 0)
        assert(max(blacklist_pixels) < nrow)
    # Find nonzero values
    row_nnz_indices, col_nnz_indices = coo.nonzero()
    row_nnz = np.bincount(row_nnz_indices, minlength=nrow)
    col_nnz = np.bincount(col_nnz_indices, minlength=ncol)
    assert(np.all(row_nnz == col_nnz))
    # Find sums
    row_sums = np.array(coo.sum(0))[0, :]
    col_sums = np.array(coo.sum(1))[:, 0]
    assert(np.all(row_sums == col_sums))
    # Return data
    return(row_nnz, row_sums)


def get_mad_filter(row_sums, mad_max):
    # Calculate MAD for log of non-zero values
    log_counts = np.log(row_sums)
    log_nnz = log_counts[row_sums > 0]
    median_log_nnz = np.median(log_nnz)
    mad_log = np.median(np.abs(log_nnz - median_log_nnz))
    # Calculate limit and generate filter
    mad_limit = median_log_nnz - (mad_log * mad_max)
    mad_filter = log_counts >= mad_limit
    return mad_filter


def get_clr_filter(clr, min_nnz, min_count, mad_max):
    # Get non-zero counts and sums for rows
    row_nnz, row_sums = get_clr_counts(clr)
    # Calculate non-zero filter
    nnz_filter = row_nnz >= min_nnz
    # Calculate count filter
    count_filter = row_sums >= min_count
    # Calculate mad filter
    mad_filter = get_mad_filter(row_sums, mad_max)
    # Generate final filter and return
    combined_filter = nnz_filter & count_filter & mad_filter
    return combined_filter


def get_filtered_bins(
    clr_paths, min_nnz, min_count, mad_max
):
    # Get bins and create empty filter
    clr = cooler.Cooler(clr_paths[0])
    clr_bins = clr.bins()[:][['chrom', 'start', 'end']]
    ok_filter = np.ones(clr_bins.shape[0]).astype('bool')
    # Loop through clr files
    for path in clr_paths:
        # Open cooler and check shape
        clr = cooler.Cooler(path)
        assert(clr_bins.equals(clr.bin()[:][['chrom', 'start', 'end']]))
        # Get filter
        clr_filter = get_clr_filter(
            clr, min_nnz=min_nnz, min_count=min_count, mad_max=mad_max
        )
        ok_filter = ok_filter & clr_filter
    # Return filter
    filtered_bins = clr_bins.loc[~ok_filter]
    return filtered_bins


# Create argument parser
parser = argparse.ArgumentParser(
    'create BED file of sparse bins'
)
parser.add_argument(
    "--clrs", nargs='+', required=True,
    help="paths to input .cool files"
)
parser.add_argument(
    "--blacklist", required=True,
    help="output bed file containing blacklisted regions"
)
parser.add_argument(
    "--min-nnz", required=True, type=int,
    help="minimum unique interactions"
)
parser.add_argument(
    "--min-count", required=True, type=int,
    help="minimum total interactions"
)
parser.add_argument(
    "--mad-max", required=True, type=float,
    help="multiples of the median absolute deviation"
)
args = parser.parse_args()
# Generate and save bins
filtered_bins = get_filtered_bins(
    args.clrs, min_nnz=args.min_nnz, min_count=args.min_count,
    mad_max=args.mad_max
)
filtered_bins.to_csv(
    args.bed, sep='\t', header=False, index=False
)
