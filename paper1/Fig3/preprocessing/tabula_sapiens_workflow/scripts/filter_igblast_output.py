import sys

import argparse
import pandas as pd
import numpy as np

############################   PARSER   ###############################

parser = argparse.ArgumentParser()
parser.add_argument('input_path', help="path to airr file")
parser.add_argument('-outdir', default='.')
parser.add_argument('--verbose', required=False, default=False,
                    dest='verbose', action='store_true',
                    help='verbose output')
parser.add_argument('-max_log_v_evalue', metavar='float_value', type=float, required=False,
                    default=-60)
parser.add_argument('-max_log_j_evalue', metavar='float_value', type=float, required=False,
                    default=-10)
parser.add_argument('-min_v_sequence_length', metavar='float_value', type=float, required=False,
                    default=160)
parser.add_argument('-min_j_sequence_length', metavar='float_value', type=float, required=False,
                    default=20)
parser.add_argument('--allow_unproductive', required=False,
                    default=False, action = 'store_true')
parser.add_argument('--allow_missing_cdr3', required=False,
                    default=False, action = 'store_true')
parser.add_argument('--allow_ns_in_sequence', required=False,
                    default=False, action = 'store_true')

args = parser.parse_args()

filename = args.input_path
outdir = args.outdir
verbose=args.verbose

MAX_LOG_V_EVALUE=args.max_log_v_evalue
MAX_LOG_J_EVALUE=args.max_log_j_evalue
ALLOW_UNPRODUCTIVE=args.allow_unproductive
ALLOW_MISSING_CDR3=args.allow_missing_cdr3
ALLOW_Ns_IN_SEQUENCE=args.allow_ns_in_sequence
MIN_V_SEQUENCE_LENGTH=args.min_v_sequence_length
MIN_J_SEQUENCE_LENGTH=args.min_j_sequence_length

sample = filename.split("/")[-1].split(".tsv.gz")[0]

###########################################################################

df = pd.read_table(filename)



vdj_pass_out_filename="{path}/{sample}_vdj.tsv.gz".format(path=outdir,sample=sample)

failed_vdj_align_out=open("{path}/{sample}_failed_vdj_alignment.fasta".format(path=outdir,sample=sample), 'w')
unproductive_out=open("{path}/{sample}_unproductive.fasta".format(path=outdir,sample=sample), 'w')
ambiguous_bases_out=open("{path}/{sample}_ambiguous_bases.fasta".format(path=outdir,sample=sample), 'w')

if __name__ == '__main__':
    
    # verify that dataframe contains anticipated columns
    for col in ['sequence',
                'v_support',
                'j_support',
                'productive',
                'cdr3',
                'v_sequence_start',
                'v_sequence_end',
                'j_sequence_start',
                'j_sequence_end']:
        if not(col in df.columns):
            raise KeyError('The following columns were not found'
                           'in dataframe: {col}'.format(col=col))
    if verbose:
        print(" Input dataframe contains {n} sequences".format(n=df.shape[0]))

    good_v_map = np.log(df.v_support.astype(float)) < MAX_LOG_V_EVALUE
    good_j_map = np.log(df.j_support.astype(float)) < MAX_LOG_J_EVALUE
    good_v_len = (df.v_sequence_end - df.v_sequence_start + 1) > MIN_V_SEQUENCE_LENGTH
    good_j_len = (df.j_sequence_end - df.j_sequence_start + 1) > MIN_J_SEQUENCE_LENGTH

    if verbose:
        if not good_v_map.all():
            print("   {n} sequences discarded because of poor V gene support.".format(
                   n=((~good_v_map).sum())))
        else:
            print("   all sequences have good V gene support.")

        if not good_j_map.all():
            print("   {n} sequences discarded because of poor J gene support.".format(
                    n=((~good_j_map).sum())))
        else:
            print("   all sequences have good J gene support.")

        if not good_v_len.all():
            print("   {n} sequences discarded because their V gene alignment is too short.".format(
                    n=((~good_v_len).sum())))
        else:
            print("   all sequences have a long enough V gene alignment.")

        if not good_j_len.all():
            print("   {n} sequences discarded because their J gene alignment is too short.".format(
                    n=((~good_j_len).sum())))
        else:
            print("   all sequences have a long enough J gene alignment.")

    good_vj_alignment = good_v_map & good_j_map & good_j_len & good_j_len

    failed_vdj_align = df[~good_vj_alignment]
    failed_vdj_align_fasta_records = ">" + failed_vdj_align.sequence_id \
                                + "|SAMPLE={}".format(failed_vdj_align.sample_id) \
                                + "\n" + failed_vdj_align.sequence + "\n"
    
    for record in failed_vdj_align_fasta_records.values:
        failed_vdj_align_out.write(record)
    failed_vdj_align_out.close()

    df = df[good_vj_alignment]


    if ALLOW_UNPRODUCTIVE:
        pass
    else:
        productive = df.productive == 'T'
        if verbose:
            if not productive.all():
                print("   {n} sequences discarded because they are not productive.".format(
                      n=((~productive).sum())))
            else:
                print("   all sequences appear productive.")

        unproductive_sequences = df[~productive]
        unproductive_sequences_fasta_records = ">" + unproductive_sequences.sequence_id \
                                    +"|SAMPLE={}".format(unproductive_sequences.sample_id) \
                                    + "\n" + unproductive_sequences.sequence + "\n"

        for record in unproductive_sequences_fasta_records.values:
            unproductive_out.write(record)

        df = df[productive]

    if ALLOW_MISSING_CDR3:
        pass
    else:
        has_cdr3 = df.cdr3.notna()

        if verbose:
            if not has_cdr3.all():
                print("   {n} sequences discarded because they do not contain a CDR3.".format(
                    n=((~has_cdr3).sum())))
            else:
                print("   all sequences appear to have a CDR3.")

        no_cdr3 = df[~has_cdr3]
        no_cdr3_fasta_records = ">" + no_cdr3.sequence_id \
                                    + "|SAMPLE={}".format(no_cdr3.sample_id) \
                                    + "\n" + no_cdr3.sequence + "\n"

        for record in no_cdr3_fasta_records.values:
            unproductive_out.write(record)
        unproductive_out.close()

        df = df[has_cdr3]

    if ALLOW_Ns_IN_SEQUENCE:
        pass
    else:
        N_in_sequence = df.sequence.map(lambda x: "N" in x)

        if verbose:
            if N_in_sequence.any():
                print("   {n} sequences discarded because they contain an N.".format(
                     n=((N_in_sequence).sum())))
            else:
                print("   all sequences have unambiguous bases.")

        ambiguous_seqs = df[N_in_sequence]
        ambiguous_seqs_fasta_records = ">" + ambiguous_seqs.sequence_id \
                                    + "|SAMPLE={}".format(ambiguous_seqs.sample_id) \
                                    + "\n" + ambiguous_seqs.sequence + "\n"

        for record in ambiguous_seqs_fasta_records.values:
            ambiguous_bases_out.write(record)
        ambiguous_bases_out.close()

        df = df[~N_in_sequence]

    if verbose:
        print(" Filtered dataframe contains {n} sequences.".format(
                n=df.shape[0]))
    
    if verbose:
        print( "Parsing sequence_id for filtered records...")
"""
    #Parse sequence_id field and drop column, seems that each sequence_id has a particular format that I'm not adhereing to
    sequence_id_records=df.sequence_id.str.split("|").map(lambda x: {item.split("=")[0]: item.split("=")[1] 
                                                                for item in x})
    sequence_id_fields = list(sequence_id_records.iloc[0].keys())
    
    rename_fields={"BARCODE":'umi_seqs',
                   "CONSCOUNT":'reads_per_umi',
                   "DUPCOUNT":'total_umis'}

    for field in sequence_id_fields:
        df[rename_fields.get(field,field)] = sequence_id_records.map(lambda x: x[field])
    
    # if all fields passing drop seq id column
    all_seq_id_fields_passing = True

    try:
        df['total_umis'] = df['total_umis'].astype(int)
    except KeyError:
        all_seq_id_fields_passing = False
        sys.stderr.write("Warning: 'information about total umis not found."
                        "(Expected field named DUPCOUNT in sequence_id)\n")
        df['total_umis'] = numpy.nan

    try:
        df['total_reads'] = df['reads_per_umi'].map(lambda x: sum(int(y) for y in x.split(",")))
    except KeyError:
        all_seq_id_fields_passing = False
        sys.stderr.write("Warning: 'information about total reads per umi not found."
                        "(Expected field named CONSCOUNT in sequence_id)\n")
        df['total_reads'] = ""

    if all_seq_id_fields_passing:
        df = df.drop('sequence_id', axis = 1)
        df['sequence_id'] = df.index

    if verbose:
        print(" Saving filtered records to {}...".format(vdj_pass_out_filename))
"""
df.to_csv(vdj_pass_out_filename, sep = '\t')

    


