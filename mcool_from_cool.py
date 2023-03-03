import argparse
import cooler
import numpy as np


def sort_cools(paths):
    '''Order cool files by bin sizes'''
    # Get bin sizes of cool files
    bin_sizes = [cooler.Cooler(path).binsize for path in paths]
    if len(bin_sizes) != len(set(bin_sizes)):
        raise ValueError('bin sizes of cool files are not unique')
    # Order paths by bin sizes and return paths and bin sizes
    path_order = np.argsort(bin_sizes)
    ordered_paths = np.array(paths)[path_order]
    ordered_sizes = np.array(bin_sizes)[path_order]
    return (ordered_paths, ordered_sizes)


def create_mcool(inpaths, outpath):
    '''Combine cool files into single mcool'''
    ordered_paths, bin_sizes = sort_cools(inpaths)
    cooler.zoomify_cooler(
        base_uris=ordered_paths, outfile=outpath, resolutions=bin_sizes,
        chunksize=1000000
    )


# Create argument parser
parser = argparse.ArgumentParser(
    'combine cool files into mcool'
)
parser.add_argument(
    '--cools', nargs='+', required=True,
    help='path to input cool files'
)
parser.add_argument(
    '--mcool', required=True,
    help='path to output cool file'
)
args = parser.parse_args()
# Create mcool
create_mcool(inpaths=args.cools, outpath=args.mcool)
