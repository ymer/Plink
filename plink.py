from __future__ import division
import bisect
import pandas as pd


def read_bed(fname, indco, snps):

    byte_translator = {}
    for byte in range(0, 256):
        binary_list = tuple((byte >> i) & 1 for i in range(8))[::-1]
        pairs = [(binary_list[i*2], binary_list[i*2 + 1]) for i in range(4)][::-1]
        byte_translator[byte] = [2 - sum(pair) if pair != (0, 1) else -1 for pair in pairs]

    genos = []

    def read_bytes(number_of_bytes):
        byte = ord(f.read(1))
        for i in range(number_of_bytes):
            geno.append(byte_translator[byte][i])

    with open(fname + '.bed', 'rb') as f:
        f.seek(3)
        for _ in range(len(snps)):
            geno = []
            for block in range(indco // 4):
                read_bytes(4)
            if indco % 4:
                read_bytes(indco % 4)
            genos.append(geno)

    return genos


def read_bim(fname):
    return pd.read_csv(fname + '.bim', sep='\s+', header=None, usecols=range(4), names=['chr', 'snp', 'gen', 'pos'])


def read_fam(fname):
    return pd.read_csv(fname + '.fam', sep='\s+', header=None, names=['famID', 'indID', 'patID', 'matID', 'sex', 'pheno'])


def read_plink(fname):
    snps = read_bim(fname)
    inds = read_fam(fname)
    genos = read_bed(fname, len(inds), snps)
    return inds, snps, genos


def write_bed(fname, inds, snps, genos):

    d = {-1: '10', 0: '00', 1: '01', 2: '11'}

    with open(fname + '.bed', 'wb') as out:
        out.write(chr(int('01101100', 2)) + chr(int('00011011', 2)) + chr(int('00000001', 2)))

        for i in range(len(snps)):
            block_index = 0
            for block_index in range((len(inds) - 1) // 4):
                genostring = [genos[i][j] for j in range(block_index * 4, block_index * 4 + 4)]
                b = ''.join(d[g] for g in genostring)[::-1]
                out.write(chr(int(b, 2)))

            block_index += 1
            if block_index * 4 < len(inds):
                lack = len(inds) - block_index * 4
                genostring = [genos[i][j] for j in range(block_index * 4, block_index*4 +lack)] + [0 for j in range(4- lack)]
                b = ''.join(d[pair] for pair in genostring)[::-1]
                out.write(chr(int(b, 2)))


def write_bim(fname, snps):
    snps.to_csv(fname + '.bim', index=False, sep=' ', header=False)


def write_fam(fname, inds):
    inds.to_csv(fname + '.fam', index=False, sep=' ', header=False)


def write_plink(fname, inds, snps, geno):
    write_fam(fname, inds)
    write_bed(fname, inds, snps, geno)
    write_bim(fname, snps)


def assoc(inds, genos):
    from scipy.stats import chi2_contingency
    pvals = []
    for i, geno in enumerate(genos):
        obs = [[0, 0], [0, 0]]
        for i, g in enumerate(geno):
            if g != -1:
                obs[0][inds.pheno[i] - 1] += 2 - g
                obs[1][inds.pheno[i] - 1] += g
        pvals.append(chi2_contingency(obs, correction=False)[1])
    return pvals

