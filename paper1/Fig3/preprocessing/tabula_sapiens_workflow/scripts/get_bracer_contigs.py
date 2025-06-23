import sys
import pandas as pd

bracer_file = sys.argv[1]
outfile = sys.argv[2]

bcrdf = pd.read_csv(bracer_file, sep = '\t') 

outfile = open(outfile, 'w')

for idx, row in bcrdf.iterrows():
    # construct sequence records
    fasta_id = row['CELL']
    sequence = str(row['SEQUENCE_INPUT'])
    sample_description = row['SEQUENCE_ID']
    record = ">{}|{}\n".format(fasta_id, sample_description)
    outfile.write(record + sequence + "\n")
outfile.close()
