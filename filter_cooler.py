import argparse
import cooler
import numpy as np
import pyranges as pr


def find_bed_pixels(bed_path, clr_bins):
    '''Find pixels overlapping intervals in bed file'''
    # Generate pyranges objects and check
    bed_pr = pr.read_bed(bed_path)
    clr_bins = clr_bins.rename(
        columns={'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'}
    )
    bins_pr = pr.from_dict(clr_bins)
    if not np.all(np.isin(bed_pr.chromosomes, bins_pr.chromosomes)):
        raise ValueError('Unknown chromosomes in blacklist')
    # Find overlapping pixels and return
    overlaps = bins_pr.count_overlaps(bed_pr)
    bed_pixels = np.where(overlaps.NumberOverlaps)[0]
    return bed_pixels


def filter_pixels(chunk, filter):
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
        chunk = filter_pixels(chunk, diagonal_filter)
    return chunk


def filter_cis(chunk, cis_only):
    '''Filter chunk pixels by removing cis interactions'''
    # Process data if pixels are to be removed
    if cis_only:
        # Extract pixel chromosomes
        bin_chroms = np.array(chunk['bins']['chrom']).astype('str')
        chrom1 = bin_chroms[chunk['pixels']['bin1_id']]
        chrom2 = bin_chroms[chunk['pixels']['bin2_id']]
        cis_filter = (chrom1 == chrom2)
        # Apply filter to chunk
        removing = (~cis_filter).sum()
        print('removing {} pixels by cis'.format(removing))
        chunk = filter_pixels(chunk, cis_filter)
    return chunk


def filter_chroms(chunk, chroms):
    '''Filter chunk pixels by removing unwanted chromosomes'''
    # Process data if pixels are to be removed
    if chroms is not None:
        # Create chromosome filter
        bin_chroms = np.array(chunk['bins']['chrom']).astype('str')
        chrom1 = bin_chroms[chunk['pixels']['bin1_id']]
        chrom2 = bin_chroms[chunk['pixels']['bin2_id']]
        chrom_filter = (np.isin(chrom1, chroms) & np.isin(chrom2, chroms))
        # Apply filter to chunk
        removing = (~chrom_filter).sum()
        print('removing {} pixels by chromosome'.format(removing))
        chunk = filter_pixels(chunk, chrom_filter)
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
        chunk = filter_pixels(chunk, blacklist_filter)
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
    # Open cooler
    inclr = cooler.Cooler(inclr_path)
    # Process supplied arguments
    if first_diagonal < 0:
        raise ValueError('minimum value of first_diagonal is 0')
    if chroms is not None:
        chroms = np.array(chroms).astype('str')
        cool_chroms = np.unique(inclr.bins()[:]['chrom']).astype('str')
        if not np.all(np.isin(chroms, cool_chroms)):
            raise ValueError('unknown chromosomes specified')
    # Get blacklisted pixels
    if blacklist is not None:
        blacklist_pixels = find_bed_pixels(
            bed_path=blacklist, clr_bins=inclr.bins()[:]
        )
    else:
        blacklist_pixels = None
    # Create iterator
    pipeline = (
        cooler.parallel.split(
            inclr, include_bins=True, map=map, chunksize=1000000
        )
        .pipe(
            filter_diagonal, first_diagonal=first_diagonal
        )
        .pipe(
            filter_cis, cis_only=cis_only
        )
        .pipe(
            filter_chroms, chroms=chroms
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
        bins=inclr.bins()[:][['chrom', 'start', 'end']],
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
