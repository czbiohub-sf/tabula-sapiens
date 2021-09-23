import sys
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

tenX_df = pd.read_csv(infile, sep = '\t') 

outfile = open(outfile, 'w')

for idx, row in bcrdf.iterrows():
    # construct sequence records
    fasta_id = row['SEQUENCE_ID']
    sequence = str(row['SEQUENCE_INPUT'])
    sample_description = row['CELL']
    record = ">{}|{}\n".format(fasta_id, sample_description)
    outfile.write(record + sequence + "\n")
outfile.close()
