import argparse
import bioframe
import cooler
import cooltools
import numpy as np
import pandas as pd


def get_exp(path, name, nproc=4):
    # Open clr and generate bioframe
    clr = cooler.Cooler(path)
    view_df = bioframe.make_viewframe(clr.chromsizes)
    # Calculate expected
    cvd = cooltools.expected_cis(
        clr=clr,
        view_df=view_df,
        smooth=True,
        aggregate_smoothed=True,
        nproc=nproc
    )
    # Add distances, deduplicate and return
    cvd['separation'] = cvd['dist'] * clr.binsize
    cvd.loc[cvd['dist'] < 2, 'balanced.avg.smoothed.agg'] = np.nan
    cvd_merged = cvd.drop_duplicates(subset=['dist'])[
        ['separation', 'balanced.avg.smoothed.agg']
    ]
    cvd_merged = cvd_merged.rename(columns={'balanced.avg.smoothed.agg': name})
    return(cvd_merged)


# Create argument parser
parser = argparse.ArgumentParser(
    'plot ps curves for multiple samples'
)
parser.add_argument(
    '--clrs', required=True, nargs='+', help=(
        'input cooler files. names can be supplied via a prefix seprated '
        'from the path via a colon. e.g. name:/path/to/clr'
    )
)
parser.add_argument(
    '--counts', required=True, help=('output file for counts')
)
parser.add_argument(
    '--nproc', type=int, default=4, help=('number of processors (default=4)')
)
args = parser.parse_args()
# Get names and paths
clr_names, clr_paths = zip(*[clr.split(':') for clr in args.clrs])
if len(clr_names) != len(set(clr_names)):
    raise ValueError('repeated name')
if len(clr_paths) != len(set(clr_paths)):
    raise ValueError('repeated path')
# Get list of counts vs distance
cvd_list = [
    get_exp(path=path, name=name, nproc=args.nproc) for
    name, path in zip(clr_names, clr_paths)
]
# Merge data
merged = cvd_list[0]
for cvd in cvd_list[1:]:
    merged = pd.merge(merged, cvd, how='outer', on='separation', sort=True)
# Save file
merged.to_csv(
    args.counts, sep='\t', header=True, index=False
)
