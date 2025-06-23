import sys, os
from functools import partial

from Bio import SeqIO, bgzf
import regex

import gzip
import argparse

import time

import pandas as pd
import numpy as np

from pacbio_vdj_utils.blast_utils import *


parser = argparse.ArgumentParser(description='Wrapper for IGHC annotation')
parser.add_argument('input_file', help="airr-formatted file")
parser.add_argument('-ighc_db', required=True, help="database containing constant region sequences")
parser.add_argument('-outdir', default='.', help="output directory (default: working directory)")
parser.add_argument('-min_c_sequence_length', type=int, required=False, default=50, 
					help="default: 50")
parser.add_argument('-max_c_sequence_length', type=int, required=False, default=2000, 
					help="default: 2000")
args = parser.parse_args()

FILENAME = args.input_file
IGHC_DATABASE = args.ighc_db
outdir = args.outdir
clen_min = args.min_c_sequence_length
clen_max = args.max_c_sequence_length

samplename=FILENAME.split("/")[-1].split('_vdj.tsv.gz')[0]
tmp_file = "{}/temp_file.fasta".format(outdir)

#########################################################################################################

df = pd.read_table(FILENAME, low_memory=False)

df['sequence_id'] = df['sequence_id'].astype(str)

#use igblast call to determine the constant region start position
df['c_sequence_start'] = df.j_sequence_end.astype(int)
df['c_sequence'] = df.apply(lambda x: x.sequence[x.c_sequence_start:], axis = 1)
clen = df.c_sequence.str.len()

sys.stderr.write("Starting with %d reads...\n" % df.shape[0])

df = df[(clen <= clen_max) & (clen >= clen_min)]

sys.stderr.write("Keeping {} reads with putative IGHCseq between {} and {}nt...\n".format(df.shape[0], 
																						  clen_min, 
																						  clen_max))
additional_blastn_options = "-word_size 9 -dust no -penalty -1 -gapopen 3 -gapextend 2 "

#now pipe to blastn
chunk_size = 1000

fasta_records = ">" + df.sequence_id.astype(str) + "\n" + df.c_sequence
fasta_records = fasta_records.values
blastn_results = []

start = time.time()
sys.stderr.write("Piping IGHC sequences to blastn...\n")
for i in range(0, df.shape[0], chunk_size):
	chunk_begin = i
	chunk_end = chunk_begin + chunk_size

	if chunk_end < df.shape[0]:
		query_string = "\n".join(fasta_records[chunk_begin:chunk_end])
	else:
		query_string = "\n".join(fasta_records[chunk_begin:])

	blastn_out, blastn_err = pipe_to_blastn(query_string, IGHC_DATABASE,
											evalue="20",
											additional_blastn_options=additional_blastn_options, 
											tmp_file=tmp_file)

	blastn_results.append(return_best_match(blastn_out))
	#os.remove(tmp_file)

	if i % 10000 == 0:
		sys.stderr.write("[%ds] %d records processed...\n" % (time.time() - start, i))
blastn_results = pd.concat(blastn_results)
blastn_results.sequence_id = blastn_results.sequence_id.astype(str)
df = df.merge(blastn_results, on = "sequence_id", how = "inner")

sys.stderr.write("Aligned {n} sequences successfully.\n".format(n = df.shape[0]))

#drop allele info
df.match = df.match.str.split("*").apply(lambda x: x[0])

df = df.rename(columns = {"match"   : "c_call",
						  "pident"  : "c_pident",
						  "length"  : "c_match_length",
						  "mismatch": "c_mismatch",
						  "gapopen" : "c_gapopen",
						  "qstart"  : "c_qstart",
						  "qend"    : "c_qend",
						  "sstart"  : "c_sstart",
						  "send"    : "c_send",
						  "evalue"  : "c_evalue",
						  "btop"    : "c_btop"
						  })

df.to_csv('{}/{}_vdjc.tsv.gz'.format(outdir,samplename), sep = '\t', index = False)
