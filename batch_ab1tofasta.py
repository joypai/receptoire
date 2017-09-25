#!/usr/bin/env python3

from Bio import SeqIO
import sys
import os.path as op
import glob

data_dir, out_dir, out_fname = sys.argv[1:]
outf = op.join(out_dir, out_fname)

ab1_files = glob.glob(data_dir+"/*.ab1")
print("converting", len(ab1_files), "ab1 files to fasta")

seqs = []
for f in ab1_files:
    record = SeqIO.read(f, "abi")
    seqs.append(record)

count = SeqIO.write(seqs, outf, "fasta")

assert len(ab1_files) == count, "conversion not successful for all sequences"
print("output fasta:", op.abspath(outf))
