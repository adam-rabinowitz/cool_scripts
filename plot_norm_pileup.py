import argparse
import cooler
import cooltools
import itertools
import matplotlib.pyplot as plt
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
        cmap='vlag', square=True, cbar_kws={'label': 'log2(obs/exp)'}, **kwargs
    )
    # Set axis labels and return
    score = data['score'].unique()[0]
    ax.set_xlabel(score)
    ax.set_ylabel('')
    return(ax)


def draw_faceted_heatmap(
    log2_avgs, names, scores, limit, ncol
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
    # Add additional data
    if scores:
        plot_data['score'] = [scores[i] for i in plot_data['name']]
    else:
        plot_data['score'] = ''
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
        data=plot_data, col='name', col_wrap=ncol, height=3, aspect=1.2,
        sharex=False
    )
    fig.figure.subplots_adjust(hspace=0.4)
    # Add plots
    fig.map_dataframe(
        draw_heatmap, vmax=vmax, vmin=vmin
    )
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
        '--loops', required=True, nargs='+', help='input loop file'
    )
    parser.add_argument(
        '--plot-width', required=True, type=int, help='width of plot in pixels'
    )
    parser.add_argument(
        '--prefix', required=True, help="prefix of output files"
    )
    parser.add_argument(
        '--score-width', type=int, help='width of scored region'
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
    if args.score_width:
        if not args.score_width % 2:
            raise ValueError('score width must be an odd number')
        if (args.plot_width / 3) < args.score_width:
            raise ValueError(
                'score width must be a third, or less, of plot width'
            )
    # Parse clrs
    clr_paths, clr_names = zip(*[clr.split(':') for clr in args.clrs])
    if len(clr_names) != len(set(clr_names)):
        raise ValueError('some names are duplicated')
    # Parse loops
    loop_paths, loop_names = zip(*[loop.split(':') for loop in args.loops])
    if len(loop_names) != len(set(loop_names)):
        raise ValueError('some names are duplicated')
    # Open clrs and check resolution
    clrs = [cooler.Cooler(clr) for clr in clr_paths]
    binsize = clrs[0].binsize
    for clr in clrs[1:]:
        if clr.binsize != binsize:
            raise ValueError('incositent bin sizes')
    # Read loops and centre on bins
    loops = [pd.read_csv(lps, sep='\t') for lps in loop_paths]
    loops = [centre_loops(lps, binsize) for lps in loops]
    # Get view dataframe
    chrom_list = [lps['chrom1'].unique() for lps in loops]
    chroms = np.unique(np.concatenate(chrom_list))
    view_df = generate_view(clrs[0], chroms=chroms)
    # Calculate expected
    expected = [cooltools.expected_cis(clr, view_df=view_df) for clr in clrs]
    # Get observed / expected pileup
    flank = ((args.plot_width - 1) // 2) * binsize
    stacks = [
        cooltools.pileup(
            clr, features_df=lps, view_df=view_df, expected_df=exp, flank=flank
        ) for
        (clr, exp), lps in itertools.product(
            zip(clrs, expected), loops
        )
    ]
    # Generate names
    names = ['\n'.join(x) for x in itertools.product(clr_names, loop_names)]
    names = [n.replace('_', ' ') for n in names]
    # Calculate log2 average of observed vs expected
    log2_avgs = [calc_log2_avg(stack) for stack in stacks]
    # Calculate scores and create plot titles
    if args.score_width:
        scores = [calc_corner_score(la, args.score_width) for la in log2_avgs]
        scores = [
            'score={:.2f} for px={}'.format(
                score, args.score_width
            ) for score in scores
        ]
    else:
        scores = None
    # Create figure and save
    fig = draw_faceted_heatmap(
        log2_avgs, names=names, scores=scores, limit=args.limit, ncol=args.ncol
    )
    fig.savefig(args.prefix + '.png', dpi=300)

