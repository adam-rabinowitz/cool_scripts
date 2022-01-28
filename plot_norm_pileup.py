import argparse
import cooler
import cooltools
import numpy as np
import pandas as pd
import seaborn as sns


def generate_view(
    clr, chroms=None
):
    # Process chromosomes
    if chroms is None:
        chroms = clr.chromnames
    else:
        for chrom in chroms:
            if chrom not in clr.chromnames:
                raise ValueError('unrecognised chromosome {}'.format(chrom))
    # Create view and return
    regions = [
        (chrom, 0, length) for
        chrom, length in zip(clr.chromnames, clr.chromsizes) if
        chrom in chroms
    ]
    view_df = pd.DataFrame(regions, columns=['chrom', 'start', 'end'])
    view_df['name'] = view_df['chrom']
    return(view_df)


def get_exp(path, chroms=None):
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
    view_df = pd.DataFrame(regions, columns=['chrom', 'start', 'end'])
    view_df['name'] = view_df['chrom']
    # Get contact vs distance and remove unwanted diagonals
    return(exp)


def centre_loops(
    loops, binsize
):
    # Get middle of bins
    middle1 = (loops['start1'] + loops['end1']) / 2
    middle2 = (loops['start2'] + loops['end2']) / 2
    # Assign middle to bins
    bin1 = np.round((middle1 - (0.5 * binsize)) / binsize).astype(np.uint32)
    bin2 = np.round((middle2 - (0.5 * binsize)) / binsize).astype(np.uint32)
    # Reset start and end columns to bins
    loops['start1'] = bin1 * binsize
    loops['start2'] = bin2 * binsize
    loops['end1'] = loops['start1'] + binsize
    loops['end2'] = loops['start2'] + binsize
    # Format and return
    loops = loops[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]
    return(loops)


def calc_log2_avg(
    stack
):
    # Calculate average
    avg = np.nanmean(stack, axis=2)
    log2_avg = np.log2(avg)
    return(log2_avg)


def calc_corner_score(
    x, pixels
):
    # Check arguments
    nrow, ncol = x.shape
    if nrow != ncol:
        raise ValueError('unexpected array shape')
    if (nrow % 2) != 1:
        raise ValueError('unexpected array size')
    if (pixels % 2) != 1:
        raise ValueError('unexpected pixel size')
    # Get centre mean
    slice_start = int((nrow - pixels) / 2)
    center_slice = slice(slice_start, slice_start + pixels)
    center = x[center_slice, center_slice]
    # Get corner mean calculate score and return
    corner = x[:pixels, -pixels:]
    score = np.nanmean(center) - np.nanmean(corner)
    return(score)


def draw_heatmap(**kwargs):
    # Extract plot datas
    data = kwargs.pop('data')
    plot_data = data.pivot(index='row', columns='column', values='value')
    # Create plot and return
    ax = sns.heatmap(
        plot_data, xticklabels=False, yticklabels=False, cbar=True,
        cmap='vlag', square=True, **kwargs
    )
    ax.set_ylabel('')    
    ax.set_xlabel('')
    return(ax)


def draw_faceted_heatmap(
    log2_avgs, names, limit, ncol
):
    # Generate plot data
    log2_avgs = np.array(log2_avgs)
    index = pd.MultiIndex.from_product(
        [range(s) for s in log2_avgs.shape], names=['name', 'row', 'column']
    )
    plot_data = pd.DataFrame(
        {'value': log2_avgs.flatten()}, index=index
    )['value']
    plot_data = plot_data.reset_index()
    plot_data['name'] = [names[i] for i in plot_data['name']]
    # Set limits
    if limit is None:
        abs_values = np.absolute([plot_data['value']])
        vmax = np.max(abs_values[np.isfinite(abs_values)])
    else:
        vmax = np.absolute(limit)
    vmin = 0 - vmax
    # Create facet grid
    fig = sns.FacetGrid(
        data=plot_data, col='name', col_wrap=ncol, height=3, aspect=1
    )
    # Add plots
    fig.map_dataframe(
        draw_heatmap, vmax=vmax, vmin=vmin
    )
    # Rename titles and return
    fig.set_titles('{col_name}')
    return(fig)


if __name__=="__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        'Plot contact enrichment for loops'
    )
    parser.add_argument(
        '--clrs', required=True, nargs='+', help='input cooler file'
    )
    parser.add_argument(
        '--loops', required=True, help='input loop file'
    )
    parser.add_argument(
        '--plot-width', required=True, type=int, help='width of plot in pixels'
    )
    parser.add_argument(
        '--score-width', required=True, type=int, help='width of scored region'
    )
    parser.add_argument(
        '--prefix', required=True, help="prefix of output files"
    )
    parser.add_argument(
        '--ncol', type=int, default=2, help="number of columns in plot"
    )
    parser.add_argument(
        '--limit', type=float, default=None, help="limit of colur scale"
    )
    args = parser.parse_args()
    # Check arguments
    if not args.plot_width % 2:
        raise ValueError('plot width must be an odd number')
    if not args.score_width % 2:
        raise ValueError('score width must be an odd number')
    if (args.plot_width / 3) < args.score_width:
        raise ValueError('score width must be a third or less of plot width')
    # Parse clrs
    clr_paths = [clr.split(':')[0] for clr in args.clrs]
    clr_names = [clr.split(':')[1] for clr in args.clrs]
    if len(clr_names) != len(set(clr_names)):
        raise ValueError('some names are duplicated')
    # Open clrs and check resolution
    clrs = [cooler.Cooler(clr) for clr in clr_paths]
    binsize = clrs[0].binsize
    for clr in clrs[1:]:
        if clr.binsize != binsize:
            raise ValueError('incositent bin sizes')
    # Read loops and centre on bins
    loops = pd.read_csv(args.loops, sep='\t')
    loops = centre_loops(loops, binsize)
    # Get view dataframe
    chroms = loops['chrom1'].unique()
    view_df = generate_view(clr, chroms=chroms)
    # Calculate expected
    expected = [cooltools.expected_cis(clr, view_df=view_df) for clr in clrs]
    # Get observed / expected pileup
    flank = ((args.plot_width - 1) // 2) * binsize
    stacks = [
        cooltools.pileup(
            clr, features_df=loops, view_df=view_df, expected_df=exp,
            flank=flank
        ) for
        clr, exp in zip(clrs, expected)
    ]
    # Calculate log2 average of observed vs expected
    log2_avgs = [calc_log2_avg(stack) for stack in stacks]
    # Calculate scores and create plot titles
    scores = [calc_corner_score(la, args.score_width) for la in log2_avgs]
    names = [
        '{}\nscore={:.2f} for px={}'.format(name, score, args.score_width) for
        name, score in zip(clr_names, scores)
    ]
    # Create figure and save
    fig = draw_faceted_heatmap(
        log2_avgs, names=names, limit=args.limit, ncol=args.ncol
    )
    fig.savefig(args.prefix + '.png', dpi=300)
