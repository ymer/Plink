# Plink
Manipulate plink files from python

Reads plink files to python. Can read .ped files, binary files and dosage files.

Example:

inds, snps = plink.read_plink('filename')

print len(snp for snp in snps if snp.chr == 1)