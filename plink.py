from __future__ import division
import bisect
import gzip
from glob import glob


class Snp:
    def __init__(self, rsNumber='?', chr='-1', pos=-1, allel1='?', allel2='?'):
        self.chr = chr
        self.rsNumber = rsNumber
        self.pos = int(pos)
        self.allel1 = allel1
        self.allel2 = allel2
        self.geno = []
        self.dosage = []

    def __repr__(self):
        return self.rsNumber

    def __cmp__(self, other):
        return position_compare(self, other)


class Ind:
    def __init__(self, famID, indID, aff, sex='0'):
        self.famID = famID
        self.indID = indID
        self.ID = famID + '_' + indID
        self.aff = aff
        self.sex = sex

    def __repr__(self):
        return self.indID


def get_byte_translator():

    byte_translator = {}
    for byte in range(0, 256):
        binary_list = tuple((byte >> i) & 1 for i in range(8))[::-1]
        pairs = [(binary_list[i*2], binary_list[i*2 + 1]) for i in range(4)][::-1]
        byte_translator[byte] = [2 - sum(pair) if pair != (0, 1) else -1 for pair in pairs]

    return byte_translator


def read_bed(fname, indco, snps, byte_translator):

    def read_bytes(number_of_bytes):
        byte = ord(f.read(1))
        for i in range(number_of_bytes):
            snp.geno.append(byte_translator[byte][i])

    with open(fname + '.bed', 'rb') as f:
        f.seek(3)
        for snp in snps:
            for block in range(indco // 4):
                read_bytes(4)
            if indco % 4:
                read_bytes(indco % 4)

    return snps


def read_bim(fname):
    return [Snp(line[1], line[0], int(line[3]), line[4], line[5]) for line in read(fname + '.bim')]


def read_fam(fname):
    return [Ind(line[0], line[1], int(line[5]), line[4]) for line in read(fname + '.fam')]


def read_plink(fname):
    byte_translator = get_byte_translator()
    snps = read_bim(fname)
    inds = read_fam(fname)
    snps = read_bed(fname, len(inds), snps, byte_translator)
    return inds, snps


def write_bed(fname, inds, snps):

    d = {-1: '10', 0: '00', 1: '01', 2: '11'}

    out = open(fname + '.bed', 'wb')
    out.write(chr(int('01101100', 2)) + chr(int('00011011', 2)) + chr(int('00000001', 2)))

    for snp in snps:
        block_index = 0
        for block_index in range((len(inds) - 1) // 4):
            genostring = [snp.geno[i] for i in range(block_index * 4, block_index * 4 + 4)]
            b = ''.join(d[g] for g in genostring)[::-1]
            out.write(chr(int(b, 2)))

        block_index += 1
        if block_index * 4 < len(inds):
            lack = len(inds) - block_index * 4
            genostring = [snp.geno[i] for i in range(block_index * 4, block_index*4 +lack)] + [0 for i in range(4- lack)]
            b = ''.join(d[pair] for pair in genostring)[::-1]
            out.write(chr(int(b, 2)))

    out.close()


def write_bim(fname, snps):
    write(fname + '.bim', [[str(snp.chr), snp.rsNumber, '-1', str(snp.pos), snp.allel1, snp.allel2] for snp in snps])


def write_fam(fname, inds):
    write(fname + '.fam', [[ind.famID, ind.indID, '0', '0', str(ind.sex), str(ind.aff)] for ind in inds])


def write_map(fname, snps):
    write(fname + '.map', [[snp.chr, snp.rsNumber, '-1', snp.pos] for snp in snps])


def write_plink(fname, inds, snps):
    write_fam(fname, inds)
    write_bed(fname, inds, snps)
    write_bim(fname, snps)


def write_dosage(fn, inds, snps):

    with open(fn + '.dosage', 'w') as out:  # Should use write function instead
        out.write('SNP\tA1\tA2\t' + '\t'.join(ind.famID + '\t' + ind.indID for ind in inds) + '\n')
        for snp in snps:
            out.write(snp.rsNumber + '\t' + snp.allel1 + '\t' + snp.allel2 + '\t' + '\t'.join(str(snp.dosage[i][0]) + '\t' + str(snp.dosage[i][1]) for i in range(len(snp.dosage))) + '\n')

    write_fam(fn + '.dosage', inds)
    write_map(fn + '.dosage', snps)


def get_edge_snps(region, snps):

    left_snp = bisect.bisect_right(snps, Snp('', region.chr, region.start - 1))
    right_snp = bisect.bisect_left(snps, Snp('', region.chr, region.end + 1)) - 1
    if left_snp <= right_snp:
        return left_snp, right_snp


def read_bim_region(fname, start, end):

    return [Snp(line[1], line[0], int(line[3]), line[4], line[5]) for i, line in enumerate(read(fname + '.bim'))
            if start <= i <= end]


def read_bed_region(fname, indco, snps, startindex, endindex, byte_translator):

    def read_bytes(number_of_bytes):
        byte = ord(f.read(1))
        for i in range(number_of_bytes):
            snp.geno.append(byte_translator[byte][i])

    if indco % 4 == 0:
        extra = 0
    else:
        extra = 1

    f = open(fname + '.bed', 'rb')
    try:
        f.seek(3 + startindex * (indco // 4 + extra))
        for snp_index, snp in enumerate(snps):
            if snp_index <= endindex:
                for block in range(indco // 4):
                    read_bytes(4)
                if indco % 4:
                    read_bytes(indco % 4)
    finally:
        f.close()

    return snps


def read_plink_region_index(fname, startsnp, endsnp):

    snps = read_bim_region(fname, startsnp, endsnp)
    inds = read_fam(fname)
    byte_translator = get_byte_translator()
    snps = read_bed_region(fname, len(inds), snps, startsnp, endsnp, byte_translator)

    return inds, snps


def read_plink_region(fname, region):

    snps = read_bim(fname)
    edges = get_edge_snps(region, snps)
    if edges:
        return read_plink_region_index(fname, *edges)
    else:
        return None, None


def position_compare(self, other):
    if self.chr < other.chr:
        return -1
    elif self.chr > other.chr:
        return 1
    else:
        if self.pos < other.pos:
            return -1
        elif self.pos > other.pos:
            return 1
        else:
            return 0


def read_map(fn):
    return [Snp(rsNumber=line[1], chr=line[0], pos=int(line[3])) for line in read(fn + '.map')]


def write(path, lines, first_line=None):

    with open(path, 'w') as f:
        if first_line:
            f.write('\t'.join(first_line) + '\n')
        f.writelines('\t'.join(str(item) for item in line) + '\n' for line in lines)


def read(path, skip_first=False, sep=''):
    with open(path) as f:
        if sep:
            rlist = [line.strip().split(sep) for line in f]
        else:
            rlist = [line.strip().split() for line in f]
        if skip_first:
            return rlist[1:]
        else:
            return rlist


def read_gzip(fn, skip_first=False, sep=False):
    with gzip.open(fn + '.gz') as f:
        if sep:
            rlist = [line.strip().split(sep) for line in f]
        else:
            rlist = [line.strip().split() for line in f]
        if skip_first:
            return rlist[1:]
        else:
            return rlist


def read_dosage(fn, gzipped=False):

    inds = read_fam(fn)
    if not gzipped:
        rows = read(fn)
    else:
        rows = read_gzip(fn)
    snps = {}
    for row in rows[1:]:
        snp = Snp(row[0], allel1=row[1], allel2=row[2])
        for i in range(3, len(row), 2):
            snp.dosage.append([float(row[i]), float(row[i+1])])
            snp.geno.append(snp.dosage[-1][0] * 2 + snp.dosage[-1][1])
        snps[snp.rsNumber] = snp

    map_snps = read_map(fn)
    for snp in map_snps:
        snps[snp.rsNumber].pos = snp.pos
        snps[snp.rsNumber].chr = snp.chr

    snps = [snp for snp in snps.values()]
    snps.sort()

    return inds, snps


def read_dosage_snps(fn, use_rs, gzipped=True):

    fns = [fn[:-4] for fn in glob(fn + '*' + '.map')]
    inds = read_fam(fns[0])

    use_rs = set(use_rs)
    snps = {}
    for i, fn in enumerate(fns):
        uselines = set()
        map_snps = read_map(fn)
        for i, snp in enumerate(map_snps):
            if snp.rsNumber in use_rs:
                uselines.add(i)
        if uselines:
            try:
                if gzipped:
                    inp = gzip.open(fn + '.gz')
                else:
                    inp = open(fn)
                rows = inp.readlines()
                rows[0] = rows[0].split()
                for index in uselines:
                    row = rows[index+1].split()
                    if row[0].upper() in use_rs:
                        snp = Snp(row[0], allel1=row[1], allel2=row[2])
                        for i in range(3, len(row), 2):
                            snp.dosage.append([float(row[i]), float(row[i+1])])
                            snp.geno.append(snp.dosage[-1][0] * 2 + snp.dosage[-1][1])
                        snps[snp.rsNumber] = snp
            finally:
                inp.close()

            for snp in map_snps:
                if snp.rsNumber in use_rs:
                    snps[snp.rsNumber].pos = snp.pos
                    snps[snp.rsNumber].chr = snp.chr

    return inds, snps

