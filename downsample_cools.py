import argparse
import cooler
import numpy as np


def sample_chunk(chunk, frac, rng):
    '''reproducibly sample chunk of cool file'''
    if frac < 1:
        # Reproducibly sample pixels
        chunk['pixels']['count'] = rng.binomial(
            chunk['pixels']['count'], frac
        )
        # Remove zero counts
        pos_count = chunk['pixels']['count'] > 0
        chunk['pixels'] = {
            key: chunk['pixels'][key][pos_count] for
            key in chunk['pixels']
        }
    return chunk


def sample_clr(inclr_path, outclr_path, count, seed):
    '''reproducibly sample cool file'''
    # Calculate fraction
    inclr = cooler.Cooler(inclr_path)
    frac = count / inclr.info['sum']
    # Create random number generator
    rng = np.random.Generator(np.random.PCG64(seed=seed))
    # Create iterator
    pipeline = (
        cooler.parallel.split(
            inclr, include_bins=False, map=map, chunksize=1000000
        )
        .pipe(
            sample_chunk, frac=frac, rng=rng
        )
        .pipe(
            lambda x: x['pixels']
        )
    )
    # Filter cooler
    cooler.create_cooler(
        outclr_path,
        bins=inclr.bins()[:],
        pixels=iter(pipeline),
        ordered=True
    )


# Create argument parser
parser = argparse.ArgumentParser(
    'reproducibly downsample cool files to their minimum count'
)
parser.add_argument(
    '--inclrs', nargs='+', required=True,
    help="input cooler files"
)
parser.add_argument(
    '--outclrs', nargs='+', required=True,
    help="output cooler files"
)
parser.add_argument(
    '--seed', type=int, default=42,
    help='seed for reproducible sampling'
)
args = parser.parse_args()
# Check arguments
if len(args.inclrs) != len(args.outclrs):
    raise ValueError('must be same number of input and output cooler files')
if len(args.inclrs) != len(set(args.inclrs)):
    raise ValueError('input cool files must be unique')
if len(args.outclrs) != len(set(args.outclrs)):
    raise ValueError('output cool files must be unique')
# Get minimum count
clr_counts = [cooler.Cooler(clr).info['sum'] for clr in args.inclrs]
min_count = min(clr_counts)
print('minimum count: {}'.format(min_count))
# Sample counts
for inclr, outclr in zip(args.inclrs, args.outclrs):
    sample_clr(
        inclr_path=inclr, outclr_path=outclr, count=min_count, seed=args.seed
    )
