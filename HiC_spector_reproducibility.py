# This file is based on the one in the HiC-spector github repository
# https://github.com/gersteinlab/HiC-spector/blob/master/run_reproducibility_v2.py
# I have altered it to enable reading directly from cool files, and allow
# multiple comparisons and tidied up a bit. If results are used in a paper
# then the following paper should be referenced
# https://doi.org/10.1093/bioinformatics/btx152
import argparse
import cooler
import itertools
import numpy as np
import pandas as pd
from scipy.sparse.linalg import eigsh


def read_matrix(cool_path, chrom):
    clr = cooler.Cooler(cool_path)
    lil_matrix = clr.matrix(sparse=True, balance=False).fetch(chrom).tolil()
    return lil_matrix


def get_Laplacian(M):
    S = M.sum(1)
    i_nz = np.where(S > 0)[0]
    S = S[i_nz]
    M = (M[i_nz].T)[i_nz].T
    S = 1 / np.sqrt(S)
    M = S*M
    M = (S*M.T).T
    n = np.size(S)
    M = np.identity(n) - M
    M = (M + M.T) / 2
    return M


def evec_distance(v1, v2):
    d1 = np.dot(v1 - v2, v1 - v2)
    d2 = np.dot(v1 + v2, v1 + v2)
    if d1 < d2:
        d = d1
    else:
        d = d2
    return np.sqrt(d)


def get_ipr(evec):
    ipr = 1.0 / (evec * evec * evec * evec).sum()
    return ipr


def get_reproducibility(M1, M2, num_evec):
    k1 = np.sign(M1.A).sum(1)
    d1 = np.diag(M1.A)
    kd1 = ~((k1 == 1) * (d1 > 0))
    k2 = np.sign(M2.A).sum(1)
    d2 = np.diag(M2.A)
    kd2 = ~((k2 == 1) * (d2 > 0))
    iz = np.nonzero((k1 + k2 > 0) * (kd1 > 0) * (kd2 > 0))[0]
    M1b = (M1[iz].A.T)[iz].T
    M2b = (M2[iz].A.T)[iz].T
    i_nz1 = np.where(M1b.sum(1) > 0)[0]
    i_nz2 = np.where(M2b.sum(1) > 0)[0]

    M1b_L = get_Laplacian(M1b)
    M2b_L = get_Laplacian(M2b)

    a1, b1 = eigsh(M1b_L, k=num_evec, which="SM")
    a2, b2 = eigsh(M2b_L, k=num_evec, which="SM")

    b1_extend = np.zeros((np.size(M1b, 0), num_evec))
    b2_extend = np.zeros((np.size(M2b, 0), num_evec))
    for i in range(num_evec):
        b1_extend[i_nz1, i] = b1[:, i]
        b2_extend[i_nz2, i] = b2[:, i]

    ipr_cut = 5
    ipr1 = np.zeros(num_evec)
    ipr2 = np.zeros(num_evec)
    for i in range(num_evec):
        ipr1[i] = get_ipr(b1_extend[:, i])
        ipr2[i] = get_ipr(b2_extend[:, i])

    b1_extend_eff = b1_extend[:, ipr1 > ipr_cut]
    b2_extend_eff = b2_extend[:, ipr2 > ipr_cut]
    num_evec_eff = min(np.size(b1_extend_eff, 1), np.size(b2_extend_eff, 1))

    evd = np.zeros(num_evec_eff)
    for i in range(num_evec_eff):
        evd[i] = evec_distance(b1_extend_eff[:, i], b2_extend_eff[:, i])

    Sd = evd.sum()
    el = np.sqrt(2)
    evs = abs(el - Sd / num_evec_eff) / el

    N = float(M1.shape[1])
    if (np.sum(ipr1 > N / 100) <= 1) | (np.sum(ipr2 > N / 100) <= 1):
        print("at least one of the maps does not look like typical Hi-C maps")
    else:
        print("size of maps: {}".format(np.size(M1, 0)))
        print("reproducibility score: {}".format(evs))
        print("num_evec_eff: {}".format(num_evec_eff))
    return evs


# Create argument parser
parser = argparse.ArgumentParser(
    'determine similarity of HiC matrices'
)
parser.add_argument(
    '--inclrs', nargs='+', required=True,
    help='input cooler files'
)
parser.add_argument(
    '--chrom', required=True,
    help='chromosome with which to determine similarity'
)
parser.add_argument(
    '--outpath', required=True,
    help='chromosome with which to determine similarity'
)
args = parser.parse_args()
# Create empty dataframe to store results
similarity = pd.DataFrame(
    data=1, columns=args.inclrs, index=args.inclrs, dtype='float'
)
# Find similarity for all combinations of cool files
for path1, path2 in itertools.combinations(args.inclrs, 2):
    print(path1)
    print(path2)
    mat1 = read_matrix(path1, args.chrom)
    mat2 = read_matrix(path2, args.chrom)
    evs = get_reproducibility(mat1, mat2, num_evec=20)
    similarity.loc[path1, path2] = evs
    similarity.loc[path2, path1] = evs
# Save similarity to filr
similarity.to_csv(
    args.outpath, sep='\t', float_format='%.4f'
)
