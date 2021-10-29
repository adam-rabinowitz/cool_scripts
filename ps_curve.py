import cooler
import cooltools.expected
import seaborn as sns


def get_cvd(clr, name, ignore_diags):
    # Specify regions if they are not supplied
    regions = [
        (chrom, 0, length) for
        chrom, length in zip(clr.chromnames, clr.chromsizes)
    ]
    # Get contact vs distance and remove unwanted diagonals
    cvd = cooltools.expected.diagsum(
        clr,
        regions=regions[0:1],
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


clr = cooler.Cooler('/g/furlong/project/90_microC/data/microC/matrices/500/mau2_0h_comb.500_cis_sampled.cool')
cvd = get_cvd(clr, 'test', 2)
def create_plot(cvd, palette):
    plot = sns.lineplot(
        data=cvd, x='dist.avg', y='balanced.avg', hue='name', palette=('red',)
    )
    plot.set(xlabel='distance', ylabel='probability')
    plot.set_xscale("log")
    plot.set_yscale("log")
    plot.legend().set_title(None)
    fig = plot.get_figure()
    fig.savefig('test.png')
