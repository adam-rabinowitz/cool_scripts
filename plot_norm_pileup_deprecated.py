import argparse
import bioframe
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
        # Circumvents bug caused by single chromosome
        if len(chroms) == 1:
            for chrom in clr.chromnames:
                if chrom not in chroms:
                    chroms = np.append(chroms, chrom)
                    break
        # Checks chromosome names
        for chrom in chroms:
            if chrom not in clr.chromnames:
                raise ValueError('unrecognised chromosome {}'.format(chrom))
    # Create view and return
    regions = [
        (chrom, 0, length, chrom) for
        chrom, length in zip(clr.chromnames, clr.chromsizes) if
        chrom in chroms
    ]
    # Create view dataframe and return
    view_df = bioframe.from_list(regions)
    return(view_df)


def read_regions(
    path
):
    # Read in data
    regions = pd.read_csv(path, sep='\t')
    # Check format and filter
    if {'chrom', 'start', 'end'}.issubset(
        regions.columns
    ):
        regions = regions[
            ['chrom', 'start', 'end']
        ]
        regions['start'] = regions['start'].astype('int64')
        regions['end'] = regions['end'].astype('int64')
    elif {'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end1'}.issubset(
        regions.columns
    ):
        regions = regions[
            ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
        ]
        regions['start1'] = regions['start1'].astype('int64')
        regions['end1'] = regions['end1'].astype('int64')
        regions['start2'] = regions['start2'].astype('int64')
        regions['end2'] = regions['end2'].astype('int64')
    else:
        raise ValueError('unexpected column names')
    # Return data
    return(regions)


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
    ax.set_xlabel('')
    ax.set_ylabel('')
    return(ax)


def draw_faceted_heatmap(
    log2_avgs, titles, limit, ncol
):
    # Generate plot data
    log2_avgs = np.array(log2_avgs)
    index = pd.MultiIndex.from_product(
        [range(s) for s in log2_avgs.shape], names=['title', 'row', 'column']
    )
    plot_data = pd.DataFrame(
        {'value': log2_avgs.flatten()}, index=index
    )['value']
    plot_data = plot_data.reset_index()
    # Add additional data
    plot_data['title'] = [titles[i] for i in plot_data['title']]
    print(plot_data)
    # Set limits
    if limit is None:
        abs_values = np.absolute([plot_data['value']])
        vmax = np.max(abs_values[np.isfinite(abs_values)])
    else:
        vmax = np.absolute(limit)
    vmin = 0 - vmax
    # Create facet grid
    plt.subplots_adjust(hspace=0.4)
    fig = sns.FacetGrid(
        data=plot_data, col='title', col_wrap=ncol, height=3, aspect=1.2
    )
    # Add plots
    fig.map_dataframe(
        draw_heatmap, vmax=vmax, vmin=vmin
    )
    fig.set_titles('{col_name}')
    return(fig)


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        'Plot contact enrichment for specified regions'
    )
    parser.add_argument(
        '--clrs', required=True, nargs='+', help='input cooler file'
    )
    parser.add_argument(
        '--regions', required=True, nargs='+', help='regions to plot'
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
        '--limit', type=float, default=None, help="limit of colour scale"
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
        raise ValueError('some clr names are duplicated')
    # Parse loops
    region_paths, region_names = zip(
        *[region.split(':') for region in args.regions]
    )
    if len(region_names) != len(set(region_names)):
        raise ValueError('some region names are duplicated')
    # Open clrs, check resolution and calculate size of flanking region
    clrs = [cooler.Cooler(clr) for clr in clr_paths]
    binsizes = [clr.binsize for clr in clrs]
    for binsize in binsizes[1:]:
        if binsize != binsizes[0]:
            raise ValueError('incositent bin sizes')
    flank = ((args.plot_width - 1) // 2) * binsizes[0]
    # Read regions and get number of columns
    regions = [read_regions(path) for path in region_paths]
    # Get columns in region files
    col_nos = [len(region.columns) for region in regions]
    for col_no in col_nos[1:]:
        if col_no != col_nos[0]:
            raise ValueError('inconsitent region file format')
    # Extract chromosomes
    if col_nos[0] == 3:
        chrom_list = [region['chrom'].unique() for region in regions]
    else:
        chrom_list = [region['chrom1'].unique() for region in regions]
    chroms = np.unique(np.concatenate(chrom_list))
    # Create view dataframe and expected counts
    view_df = generate_view(clrs[0], chroms=chroms)
    expected = [cooltools.expected_cis(clr, view_df=view_df) for clr in clrs]
    # Get observed / expected pileup and count depth
    stacks = [
        cooltools.pileup(
            clr, features_df=region, view_df=view_df, expected_df=exp,
            flank=flank
        ) for
        (clr, exp), region in itertools.product(
            zip(clrs, expected), regions
        )
    ]
    for stack in stacks:
        print(np.nanmin(stack, axis=2))
    counts = [stack.shape[2] for stack in stacks]
    # Generate names
    titles = ['\n'.join(x) for x in itertools.product(clr_names, region_names)]
    titles = [n.replace('_', ' ') for n in titles]
    # Calculate log2 average of observed vs expected
    log2_avgs = [calc_log2_avg(stack) for stack in stacks]
    # Generate metrics
    if args.score_width:
        scores = [calc_corner_score(la, args.score_width) for la in log2_avgs]
        metrics = [
            'n={} & score={:.2f} for px={}'.format(
                count, score, args.score_width
            ) for count, score in zip(counts, scores)
        ]
    else:
        metrics = ['n={}'.format(count) for count in counts]
    # Add metrics to titles
    titles = [t + '\n' + m for t, m in zip(titles, metrics)]
    # Create figure and save
    fig = draw_faceted_heatmap(
        log2_avgs, titles=titles, limit=args.limit, ncol=args.ncol
    )
    fig.savefig(args.prefix + '.png', dpi=300)
