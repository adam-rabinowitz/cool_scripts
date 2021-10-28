import argparse
import cooler
import numpy as np


def filter_clr(clr, chroms, diagonal, cis_only):
    # Get bins and pixels
    bins = clr.bins()[:]
    pixels = clr.pixels()[:]
    # Process chromosome offsets
    chrom_offsets = [clr.offset(c) for c in clr.chromnames]
    assert(chrom_offsets[0] == 0)
    assert(chrom_offsets[-1] < bins.shape[0])
    chrom_offsets = chrom_offsets[1:] + [bins.shape[0]]
    # Get chromosome indices
    if chroms is None:
        chroms = clr.chromnames
    chrom_indices = [clr.chromnames.index(c) for c in chroms]
    # Add chromosome features
    pixels['chrom1'] = np.digitize(pixels['bin1_id'], chrom_offsets)
    pixels['chrom2'] = np.digitize(pixels['bin2_id'], chrom_offsets)
    pixels['wanted_chrom1'] = pixels['chrom1'].isin(chrom_indices)
    pixels['wanted_chrom2'] = pixels['chrom2'].isin(chrom_indices)
    # Add diagonal data
    pixels['diagonal'] = pixels['bin2_id'] - pixels['bin1_id']
    assert(all(pixels['diagonal'] >= 0))
    # Filter pixels
    if cis_only:
        pixels['filter'] = (
            pixels['wanted_chrom1'] &
            (pixels['chrom1'] == pixels['chrom2']) &
            (pixels['diagonal'] >= diagonal)
        )
    else:
        pixels['filter'] = (
            pixels['wanted_chrom1'] &
            pixels['wanted_chrom2'] &
            (pixels['diagonal'] >= diagonal)
        )
    # Filter data and create cooler
    filtered_pixels = pixels.loc[
        pixels['filter'],
        ['bin1_id', 'bin2_id', 'count']
    ].reset_index(drop=True)
    return(bins, filtered_pixels)


# Create argument parser
parser = argparse.ArgumentParser(
    'calculate correlation between chromosome interaction counts'
)
parser.add_argument(
    "--incool", required=True,
    help="path to input .cool file"
)
parser.add_argument(
    "--outcool", required=True,
    help="path to output .cool file"
)
parser.add_argument(
    "--chroms", nargs='+',
    help="chromosome from which to extract counts"
)
parser.add_argument(
    "--diagonal", default=0, type=int,
    help="diagonals smaller than this number will be excluded"
)
parser.add_argument(
    "--cis-only", action='store_true',
    help="remove trans counts"
)
args = parser.parse_args()
if (args.chroms is None) & (args.diagonal == 0) & (not args.cis_only):
    parser.error('no filtering has been specified')
# Open cooler and filter
clr = cooler.Cooler(args.incool)
bins, pixels = filter_clr(
    clr, chroms=args.chroms, diagonal=args.diagonal,
    cis_only=args.cis_only
)
# Write new cooler
cooler.create_cooler(
    args.outcool, bins=bins, pixels=pixels, metadata=clr.info['metadata'],
    assembly=clr.info['genome-assembly'], symmetric_upper=True,
    ordered=True
)
