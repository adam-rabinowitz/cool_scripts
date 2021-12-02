import argparse
import cooler
import cooltools.expected
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns


def get_cvd(path, name, ignore_diags, chroms=None):
    # Open clr and process chromosomes
    clr = cooler.Cooler(path)
    if chroms is None:
        chroms = clr.chromnames
    else:
        for chrom in chroms:
            if chrom not in clr.chromnames:
                raise ValueError('unrecognised chromosome {}'.format(chrom))
    # Specify regions if they are not supplied
    regions = [
        (chrom, 0, length) for
        chrom, length in zip(clr.chromnames, clr.chromsizes) if
        chrom in chroms
    ]
    # Get contact vs distance and remove unwanted diagonals
    cvd = cooltools.expected.diagsum(
        clr,
        regions=regions,
        transforms={
            'balanced': lambda p: p['count']*p['weight1']*p['weight2']
        },
        ignore_diags=ignore_diags
    )
    # Calculate combined expected and expected
    exp, exp_slope = cooltools.expected.logbin_expected(cvd)[0:2]
    comb_exp, comb_exp_slope = cooltools.expected.combine_binned_expected(
        binned_exp=exp, binned_exp_slope=exp_slope
    )[0:2]
    # Format data
    comb_exp['name'] = name
    comb_exp['dist.avg'] = comb_exp['diag.avg'] * clr.binsize
    comb_exp = comb_exp[["name", "dist.avg", "balanced.avg"]]
    comb_exp_slope['name'] = name
    comb_exp_slope['dist.avg'] = comb_exp_slope['diag.avg'] * clr.binsize
    comb_exp_slope = comb_exp_slope[["name", "dist.avg", "slope"]]
    return(comb_exp, comb_exp_slope)


def create_exp_plot(exp, palette):
    ax = sns.lineplot(
        data=exp, x='dist.avg', y='balanced.avg', hue='name', palette=palette
    )
    ax.set(xlabel='distance', ylabel='probability')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend().set_title(None)
    fig = ax.get_figure()
    return(fig)


def create_exp_slope_plot(exp_slope, palette):
    ax = sns.lineplot(
        data=exp_slope, x='dist.avg', y='slope', hue='name', palette=palette
    )
    ax.set(xlabel='distance', ylabel='slope')
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
    '--clrs', required=True, nargs='+', help=(
        'input cooler files. names can be supplied via a prefix seprated '
        'from the path via a colon. e.g. name:/path/to/clr'
    )
)
parser.add_argument(
    '--plot', required=True, help='output plot'
)
parser.add_argument(
    '--colors', nargs='+', help='input loop file'
)
parser.add_argument(
    '--chroms', nargs='+',
    help='chromosomes over which to calculate probability'
)
parser.add_argument(
    '--ignore-diags', type=int, help='ignore diagonals below this number'
)
args = parser.parse_args()
# Get names
names = [x.split(':')[0] for x in args.clrs]
clrs = [x.split(':')[-1] for x in args.clrs]
# Check colours
if args.colors is not None:
    if len(args.colors) != len(clrs):
        raise ValueError('single color must be provided for each clr')
    colors = {x: '#' + y for x, y in zip(names, args.colors)}
else:
    colors = {x: None for x in names}
# Read data and concatenate
exp_list = [
   get_cvd(clr, name, args.ignore_diags, args.chroms) for
   clr, name in zip(clrs, names)
]
exp_df = pd.concat(
    [x[0] for x in exp_list], axis=0, ignore_index=True
)
exp_slope_df = pd.concat(
    [x[1] for x in exp_list], axis=0, ignore_index=True
)
# Create plot and save
fig, (ax1, ax2) = plt.subplots(
    nrows=2, ncols=1, sharex=True, figsize=[6.4, 7.4], dpi=200
)
fig.suptitle('P(s) Plots')
# Create P(s) curve plot
for name in exp_df['name'].unique():
    plot_data = exp_df[exp_df['name'] == name]
    ax1.plot(
        plot_data['dist.avg'], plot_data['balanced.avg'],
        label=name, color=colors[name]
    )
# Format P(s) curve plot
ax1.set_ylabel('probability')
ax1.xaxis.set_tick_params(which='both', labelbottom=True)
ax1.legend()
ax1.set_xscale("log")
ax1.set_yscale("log")
# Create slope plot
for name in exp_slope_df['name'].unique():
    plot_data = exp_slope_df[exp_slope_df['name'] == name]
    ax2.plot(
        plot_data['dist.avg'], plot_data['slope'], label=name,
        color=colors[name]
    )
# Format slope plot
ax2.set_ylabel('slope')
ax2.set_xlabel('separation')
ax2.legend()
ax2.set_xscale("log")
# Save plot
fig.savefig(args.plot, dpi=300)
