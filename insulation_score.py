import argparse
import cooler
import cooltools.insulation
import pandas as pd

def get_insulation(path, windows, ignore_diags, chroms):
    # Open clr and process chromosomes
    clr = cooler.Cooler(path)
    if chroms is not None:
        for chrom in chroms:
            if chrom not in clr.chromnames:
                raise ValueError('unrecognised chromosome {}'.format(chrom))
    # Get insulation data, flter and return
    ins = cooltools.insulation.calculate_insulation_score(
        clr, window_bp=windows, ignore_diags=ignore_diags
    )
    if chroms is not None:
        ins = ins[ins.chrom.isin(args.chroms)]
    return(ins)

# Create argument parser
parser = argparse.ArgumentParser(
    'plot ps curves for multiple samples'
)
parser.add_argument(
    '--clr', required=True, help='input cooler files'
)
parser.add_argument(
    '--windows', nargs='+', type=int, required=True,
    help='list of sliding window sizes'
)
parser.add_argument(
    '--scores', required=True, help='output file containing scores'
)
parser.add_argument(
    '--ignore-diags', type=int, default=None,
    help='ignore diagonals below this number'
)
parser.add_argument(
    '--chroms', nargs='+', default=None,
    help='chromosomes for which to report data'
)
args = parser.parse_args()
# Get insulation scores
ins = get_insulation(
    path=args.clr, windows=args.windows, ignore_diags=args.ignore_diags,
    chroms=args.chroms
)
ins.to_csv(
    args.scores, sep='\t', header=True, index=False
)
