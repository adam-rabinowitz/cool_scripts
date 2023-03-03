import argparse
import cooler
import numpy as np
import pandas as pd
import random
from scipy.stats.mstats import gmean
import seaborn as sns


def draw_heatmap(**kwargs):
    # Extract arguments
    data = kwargs.pop('data')
    index = kwargs.pop('row')
    columns = kwargs.pop('column')
    values = kwargs.pop('value')
    # Pivot data and return arguments
    plot_data = data.pivot(index=index, columns=columns, values=values)
    ax = sns.heatmap(
        plot_data, xticklabels=False, yticklabels=False, cmap='fall', **kwargs
    )
    return(ax)


def draw_faceted_heatmaps(df):
    # Create facet grid
    fig = sns.FacetGrid(
        data=df, col='name', col_wrap=4, height=3, aspect=1
    )
    # Add plots
    fig.map_dataframe(
        draw_heatmap, row='row', column='column', value='value', cbar=True,
        square=True, vmin=0
    )
    # Rename titles and return
    fig.set_titles('{col_name}')
    return(fig)


def get_slices(
    clr, bin, expand
):
    # Get extent of bin and its chromosome
    bin_start, bin_end = clr.extent(bin)
    assert(bin_end - bin_start == 1)
    bin_index = bin_start
    chrom_start, chrom_end = clr.extent(bin[0])
    # Expand bin and locate on chromosome
    expand_start = bin_start - expand
    expand_end = bin_end + expand
    possible_bins = np.arange(expand_start, expand_end, 1)
    chr_filter = (possible_bins >= chrom_start) & (possible_bins < chrom_end)
    # Determine which bins to fill in output matrix
    output_bins = np.where(chr_filter)[0]
    output_slice = slice(output_bins[0], output_bins[-1] + 1)
    # Determine which bins to extract counts for
    count_bins = possible_bins[chr_filter]
    count_slice = slice(count_bins[0], count_bins[-1] + 1)
    # Return slices
    return(bin_index, count_slice, output_slice)


def get_interaction_counts(
    clr, bin1, bin2, expand, balance=False, name=None
):
    # Get bin data
    bin1_index, count1_slice, output1_slice = get_slices(clr, bin1, expand)
    bin2_index, count2_slice, output2_slice = get_slices(clr, bin2, expand)
    # Get slice of clr matrix for requested bins
    count_array = clr.matrix(balance=balance)[count1_slice, count2_slice]
    count_array = count_array.astype(np.float32)
    # Set pixels from lower half of clr matrix to nan
    if bin1[0] == bin2[0]:
        assert(bin2_index >= bin1_index)
        offset = bin1_index - bin2_index
        lower_indices = np.tril_indices_from(count_array, k=offset)
        count_array[lower_indices] = np.nan
    # Generate output matrix and add counts
    dim = dim = (2 * expand) + 1
    output_array = np.full((dim, dim), np.nan)
    output_array[output1_slice, output2_slice] = count_array
    # Convert to dataframe
    output_df = pd.DataFrame(output_array).reset_index().melt('index')
    output_df.columns = ['row', 'column', 'value']
    # Add name and return
    if name is None:
        name = '{}:{}-{}\n{}:{}-{}'.format(*bin1, *bin2)
    output_df['name'] = name
    return(output_df)


def read_loops(
    path, columns, minsize=0, header=True, sep='\t'
):
    # Get columns
    assert(len(columns) == 6)
    bin1_cols = columns[0:3]
    bin2_cols = columns[3:6]
    # Get loop anchors
    loop_list = []
    with open(path, 'rt') as infile:
        # Skip header
        if header:
            _ = infile.readline()
        # Get remaining data
        for line in infile:
            # Extract bin data
            line_list = line.strip().split(sep)
            bin1 = [line_list[i] for i in bin1_cols]
            bin1[1], bin1[2] = int(bin1[1]), int(bin1[2])
            bin2 = [line_list[i] for i in bin2_cols]
            bin2[1], bin2[2] = int(bin2[1]), int(bin2[2])
            # filter by loop size
            loop_size = bin2[1] - bin1[2]
            assert(loop_size) > 0
            if loop_size >= minsize:
                loop_list.append([bin1, bin2])
    return(loop_list)


def nonzero_geomean(a):
    a_filt = a[~np.isnan(a) & ~(a == 0)]
    geo_mean = gmean(a_filt)
    return(geo_mean)


def calculate_average(count_list):
    # Merge data and select average
    count_df = pd.concat(count_list, axis=0, ignore_index=True)
    avg = count_df.pivot_table(
        values='value', index='row', columns='column', aggfunc=nonzero_geomean
    )
    # Reformat average and return
    avg = avg.melt(ignore_index=False)
    avg.insert(0, 'row', avg.index)
    avg.index = range(avg.shape[0])
    avg['name'] = 'average'
    return(avg)


# Create argument parser
parser = argparse.ArgumentParser(
    'reproducibly sample cool files to provided count'
)
parser.add_argument(
    'clr', help="input cooler file"
)
parser.add_argument(
    'loops', help="input loop file"
)
parser.add_argument(
    'plot', help="output plot"
)
parser.add_argument(
    '--expand', type=int, default=9,
    help='number of bins by which to expand view (default=9)'
)
parser.add_argument(
    '--minsize', type=int, default=0, help='minimum size of loops (default=0)'
)
parser.add_argument(
    '--average', action='store_true', help='plot avergae of all loops'
)
parser.add_argument(
    '--nloops', type=int, help='number of loops to plot'
)
parser.add_argument(
    '--random', action='store_true', help='randomly select loops to plot'
)
parser.add_argument(
    '--raw', action='store_true', help='plot raw counts instead of normalised'
)
parser.add_argument(
    '--columns', default='1,2,3,4,5,6',
    help='comma seperate list of columns containing chromosome, start and end for each bin'
)
args = parser.parse_args()
# Get columns
columns = [int(x) - 1 for x in args.columns.split(',')]
# Get loops
loop_list = read_loops(args.loops, columns=columns, minsize=args.minsize)
# Open cooler and get counts
clr = cooler.Cooler(args.clr)
count_list = [
    get_interaction_counts(
        clr, bin1=b1, bin2=b2, expand=args.expand, balance=not(args.raw)
    ) for
    b1, b2 in loop_list
]
# Calculate average
if args.average:
    average = calculate_average(count_list)
# Select data for plotting
if args.nloops is not None:
    # Calculate number of loops
    nloops = min(args.nloops, len(loop_list))
    # Select data for plotting
    if args.random:
        count_list = random.sample(count_list, nloops)
    else:
        count_list = count_list[:nloops]
# Merge data
if args.average:
    count_list = [average] + count_list
count_df = pd.concat(count_list, axis=0, ignore_index=True)
# Create plot
fig = draw_faceted_heatmaps(count_df)
fig.savefig(args.plot, dpi=300)

