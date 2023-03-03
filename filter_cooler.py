import argparse
import cooler
import numpy as np


def find_blacklist_pixels(bed_path, clr):
    '''Find pixels overlapping intervals in bed file'''
    # Create empty set to store pixels
    bed_pixels = set()
    # Loop through lines of bed file
    with open(bed_path) as bed:
        for line in bed:
            # Get location listed in line
            region = line.strip().split('\t')[0:3]
            pixel_range = clr.extent(region)
            bed_pixels.update(range(*pixel_range))
    # Format bed pixels and return
    blacklist_pixels = np.sort(
        np.array(list(bed_pixels)).astype('int64')
    )
    return blacklist_pixels


def apply_filter(chunk, filter):
    '''Apply filter to chunk'''
    # Apply filter to chunk and return
    chunk['pixels'] = {
        key: chunk['pixels'][key][filter] for
        key in chunk['pixels']
    }
    return chunk


def filter_diagonal(chunk, first_diagonal):
    '''Filter chunk pixels by removing unwanted diagonals'''
    # Process data if pixels are to be removed
    if first_diagonal > 0:
        # Create filter
        diagonal = chunk['pixels']['bin2_id'] - chunk['pixels']['bin1_id']
        diagonal_filter = (diagonal >= first_diagonal)
        # Apply filter to chunk
        removing = (~diagonal_filter).sum()
        print('removing {} pixels by diagonal'.format(removing))
        chunk = apply_filter(chunk, diagonal_filter)
    return chunk


def filter_cis(chunk, clr_bins, cis_only):
    '''Filter chunk pixels by removing cis interactions'''
    # Process data if pixels are to be removed
    if cis_only:
        # Extract pixel chromosomes
        bin_chroms = np.array(clr_bins['chrom']).astype('str')
        chrom1 = bin_chroms[chunk['pixels']['bin1_id']]
        chrom2 = bin_chroms[chunk['pixels']['bin2_id']]
        cis_filter = (chrom1 == chrom2)
        # Apply filter to chunk
        removing = (~cis_filter).sum()
        print('removing {} pixels by cis'.format(removing))
        chunk = apply_filter(chunk, cis_filter)
    return chunk


def filter_chroms(chunk, clr_bins, chroms):
    '''Filter chunk pixels by removing unwanted chromosomes'''
    # Process data if pixels are to be removed
    if chroms is not None:
        # Create chromosome filter
        bin_chroms = np.array(clr_bins['chrom']).astype('str')
        chrom1 = bin_chroms[chunk['pixels']['bin1_id']]
        chrom2 = bin_chroms[chunk['pixels']['bin2_id']]
        chrom_filter = (np.isin(chrom1, chroms) & np.isin(chrom2, chroms))
        # Apply filter to chunk
        removing = (~chrom_filter).sum()
        print('removing {} pixels by chromosome'.format(removing))
        chunk = apply_filter(chunk, chrom_filter)
    return chunk


def filter_blacklist(chunk, blacklist_pixels):
    '''Filter chunk pixels by removing blacklisted pixels'''
    # Create blacklist filter
    if blacklist_pixels is not None:
        # Create blacklist filter
        blacklist_filter = np.logical_and(
            np.isin(chunk['pixels']['bin1_id'], blacklist_pixels, invert=True),
            np.isin(chunk['pixels']['bin2_id'], blacklist_pixels, invert=True)
        )
        # Apply filter to chunk
        removing = (~blacklist_filter).sum()
        print('removing {} pixels by blacklist'.format(removing))
        chunk = apply_filter(chunk, blacklist_filter)
    return chunk


def get_pixels(chunk):
    '''Get pixels for output'''
    pixels = chunk['pixels']
    retaining = len(pixels['bin1_id'])
    print('retaining {} pixels'.format(retaining))
    return pixels


def filter_cool(
    inclr_path, outclr_path, first_diagonal, cis_only, chroms, blacklist
):
    '''Combine and apply individual filters'''
    # Open cooler and get bins
    inclr = cooler.Cooler(inclr_path)
    clr_bins = inclr.bins()[:]
    # Process supplied arguments
    if first_diagonal < 0:
        raise ValueError('minimum value of first_diagonal is 0')
    if chroms is not None:
        chroms = np.array(chroms).astype('str')
        cool_chroms = np.array(inclr.chromnames)
        if not np.all(np.isin(chroms, cool_chroms)):
            raise ValueError('unknown chromosomes specified')
    # Get blacklisted pixels
    if blacklist is not None:
        blacklist_pixels = find_blacklist_pixels(bed_path=blacklist, clr=inclr)
    else:
        blacklist_pixels = None
    # Create iterator
    pipeline = (
        cooler.parallel.split(
            inclr, include_bins=False, map=map, chunksize=1000000
        )
        .pipe(
            filter_diagonal, first_diagonal=first_diagonal
        )
        .pipe(
            filter_cis, clr_bins=clr_bins, cis_only=cis_only
        )
        .pipe(
            filter_chroms, clr_bins=clr_bins, chroms=chroms
        )
        .pipe(
            filter_blacklist, blacklist_pixels=blacklist_pixels
        )
        .pipe(
            get_pixels
        )
    )
    # Filter cooler
    cooler.create_cooler(
        outclr_path,
        bins=clr_bins[['chrom', 'start', 'end']],
        pixels=iter(pipeline),
        ordered=True
    )


# Create argument parser
parser = argparse.ArgumentParser(
    'filter cool files'
)
parser.add_argument(
    '--incool', required=True,
    help='path to input cool file'
)
parser.add_argument(
    '--outcool', required=True,
    help='path to output cool file'
)
parser.add_argument(
    '--chroms', nargs='+', default=None,
    help='chromosomes from which to keep pixels'
)
parser.add_argument(
    '--diagonal', default=0, type=int,
    help='pixels from smaller diagonals will be removed'
)
parser.add_argument(
    '--blacklist', default=None,
    help='blacklist bed file'
)
parser.add_argument(
    '--cis-only', action='store_true',
    help='remove trans counts'
)
args = parser.parse_args()
# Open cooler and filter
filter_cool(
    inclr_path=args.incool,
    outclr_path=args.outcool,
    first_diagonal=args.diagonal,
    cis_only=args.cis_only,
    chroms=args.chroms,
    blacklist=args.blacklist
)
