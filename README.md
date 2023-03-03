# cool scripts

This repository contains scripts for processing cool files. All scripts utilise the python package cooler (https://cooler.readthedocs.io/en/latest/index.html) and were tested with version 0.9.1 of this package.

## Script summary

This repository contains the following scripts. The arguments for all scripts can be viewed by running the script with the '-h' option.

* **downsample_cools.py** Downsamples a set of cool files so that the outputted cool files all have the same interaction count as the input cool file with the lowest count. The use of a seed for the random number generator used to sample the cool files ensures this downsampling is reproducible.

* **filter_cooler.py** Filters cool files to remove unwanted diagonals, unwanted chromosomes, interactions involving blacklisted regions and trans interactions.

## Generating directly comparable cool files

To generate cool files that are directly comparable to one another one should process the raw cool files in the following four steps. These steps generate normalised cool files with identical interactions counts across identical pixels and that have been normalised across the same bins. 

1. Filter individual cool files with the script 'filter_cooler.py'. This script allows you to remove unwanted diagonals, trans interactions, unwanted chromosomes and interactions involving blacklisted regions. To allow direct comparison between cool files each of the files should be filtered with the exact same parameters.

2. Downsample cool files using the script 'downsample_cools.py'. This script downsamples a set of input cool files so that the generated output cool files will have the same count (+/- <0.01%) as the input cool file with the lowest interaction count. 

3. Identify bins to be excluded from the normalisation process using the script 'create_balance_blacklist.py'. The output of the script is a BED file of bins that fail to pass the bin filters in one or more of the input cool files. These bins will be excluded from the normalisation process in the next step. 

4. Perform cool file normalisation using the BED file generated in the previous step. The normalisation should be perfomed with the command line function 'balance' from the cooler package. Only bins defined in the BED file should be excluded from the normalisation. This is achieved by supplying the BED file to the argument --blacklist and setting the following parameters to 0: --ignore-diags, --mad-max, --min-nnz and --min-count.

