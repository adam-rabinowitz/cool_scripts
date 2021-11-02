import argparse
import cooler
import cooltools.expected
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def get_cvd(path, name, ignore_diags):
    # Open clr
    clr = cooler.Cooler(path)
    # Specify regions if they are not supplied
    regions = [
        (chrom, 0, length) for
        chrom, length in zip(clr.chromnames, clr.chromsizes)
    ]
    # Get contact vs distance and remove unwanted diagonals
    cvd = cooltools.expected.diagsum(
        clr,
        regions=regions,
        transforms={'balanced': lambda p: p['count']*p['weight1']*p['weight2']},
        ignore_diags=ignore_diags
    )
    # Aggregate data across chromosomes and add name
    cvd_agg = cvd.groupby('diag').agg(
        {'n_valid':'sum', 'count.sum':'sum', 'balanced.sum':'sum'}
    ).reset_index()
    # Calculate binned expected and add distance
    binned_cvd = cooltools.expected.logbin_expected(cvd)[0]
    binned_cvd['name'] = name
    binned_cvd['dist.avg'] = binned_cvd['diag.avg'] * clr.binsize
    binned_cvd = binned_cvd[["name", "dist.avg", "balanced.avg"]]
    return(binned_cvd)


def create_plot(cvd, palette):
    ax = sns.lineplot(
        data=cvd, x='dist.avg', y='balanced.avg', hue='name', palette=palette
    )
    ax.set(xlabel='distance', ylabel='probability')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend().set_title(None)
    fig = ax.get_figure()
    return(fig)


# Create argument parser
parser = argparse.ArgumentParser(
    'plot ps curves for multiple samples'
)
parser.add_argument(
    '--clrs', required=True, nargs='+', help=('input cooler files. names can be '
    'supplied via a prefix seprated from the path via a colon. e.g. '
    'name:/path/to/clr')
)
parser.add_argument(
    '--plot', help='output plot'
)
parser.add_argument(
    '--colors', nargs='+', help='input loop file'
)
parser.add_argument(
    '--ignore-diags', type=int, help='ignore diagonals below this number'
)
args = parser.parse_args()
# Check arguments
if args.colors is not None:
    if len(args.colors) != len(args.clrs):
        raise ValueError('single color must be provided for each clr')
    colors = ['#' + col for col in args.colors]
else:
    colors = None
# Get names
names = [x.split(':')[0] for x in args.clrs]
clrs = [x.split(':')[-1] for x in args.clrs]
# Read data and concatenate
cvd_list = [
   get_cvd(clr, name, args.ignore_diags) for
   clr, name in zip(clrs, names)
]
cvd_df = pd.concat(cvd_list, axis=0, ignore_index=True)
# Create plot and save
cvd_plot = create_plot(cvd_df, colors)
cvd_plot.savefig(args.plot, dpi=300)

