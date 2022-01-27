import argparse
import cooler
import cooltools
import collections
import pandas as pd
import os
import pyBigWig
import re

# Function to extract chromosome lengths from bed file
def get_chromosome_lengths(
    path
):
    # Create output variables
    lengths = []
    # Loop through lines of bed file
    with open(path, 'rt') as infile:
        for line in infile:
            chrom, length = line.strip().split('\t')
            lengths.append((chrom, int(length)))
    return(lengths)

# Read insulation files
def parse_insulation_df(
    ins
):
    # Create output and get score columns
    scores = collections.OrderedDict()
    score_cols = [
        col for col in ins.columns.values if
        re.match('log2_insulation_score_\\d+', col)
    ]
    # Generate dataframe for each score column
    for col in score_cols:
        # Get resolution
        resolution = re.sub('log2_insulation_score_', '', col)
        # Select data and rename columns
        select_ins = ins[['chrom', 'start', 'end', col]]
        select_ins.columns = ['chrom', 'start', 'end', 'score']
        # Filter data and store
        select_ins = select_ins[~select_ins['score'].isna()]
        scores[resolution] = select_ins
    return(scores)

# Generates bigwig file for ligations from a single capture
def generate_bigwig(
    df, chrom_lengths, bigwig
):
    # Sort by chromosome
    chroms = [cl[0] for cl in chrom_lengths]
    df['chrom'] = pd.Categorical(df['chrom'], chroms)
    df = df.sort_values(['chrom', 'start'])
    # Open bigwig, write and close
    bw = pyBigWig.open(bigwig, 'w')
    bw.addHeader(chrom_lengths)
    bw.addEntries(
        df['chrom'].tolist(), df['start'].tolist(),
        ends=df['end'].tolist(), values=df['score'].tolist()
    )
    bw.close()

if __name__ == '__main__':
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='idenitfy ligation events from BAM/SAM file'
    )
    parser.add_argument(
        'cool', help='cool file'
    )
    parser.add_argument(
        'windows', help='comma seperate list of windows in bp'
    )
    parser.add_argument(
        'prefix', help='output file prefix'
    )
    args = parser.parse_args()
    # Generate output directory
    abs_path = os.path.abspath(args.prefix)
    abs_dir = os.path.dirname(abs_path)
    if not os.path.isdir(abs_dir):
        os.makedirs(abs_dir)
    # Open clr and get chromosome lengths
    clr = cooler.Cooler(args.cool)
    chrom_lengths = list(
        zip(clr.chroms()[:]['name'], clr.chroms()[:]['length'])
    )
    # Get insulation and save to file
    windows = list(map(int, args.windows.split(',')))
    ins = cooltools.insulation(clr, window_bp=windows)
    ins.to_csv(args.prefix + '.txt', sep='\t')
    # Parse insulation file and generate bigwigs
    scores = parse_insulation_df(ins)
    for res in scores:
        bigwig = '.'.join([args.prefix, res, 'bw'])
        generate_bigwig(scores[res], chrom_lengths, bigwig)
