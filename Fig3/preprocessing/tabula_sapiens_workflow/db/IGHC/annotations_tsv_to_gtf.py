import pandas as pd

df = pd.read_table('IGHC_annotations.tsv')

df['sequence'] = '14'
df['source'] = "IMGT"
df['feature'] = df['IMGT label'].map(lambda x: 'gene' if ('C-GENE-UNIT' in x) else 'exon')
df['start'] = df['Positions in NCBI Accession number'].map(lambda x: int(x.split("..")[0]))
df['end'] = df['Positions in NCBI Accession number'].map(lambda x: int(x.split("..")[1]))
df['score'] = '.'
df['strand'] = df["Orientation"].map(lambda x: '-' if (x == 'REV') else "+")
df['phase'] = '.'

df['attributes'] = 'gene_name='+df.Gene + df['IMGT label'].map(lambda x: '' if ('C-GENE-UNIT' in x)
                                                                             else ";exon_name="+x)

df = df[['sequence','source','feature','start','end','score','strand','phase','attributes']]
df = df[df.feature == 'gene']

df.to_csv('IGHC_genes.gff', sep='\t', index=False, header=False)
print(df)

#     sequence    The name of the sequence where the feature is located.
# 2   source  Keyword identifying the source of the feature, like a program (e.g. Augustus or RepeatMasker) or an organization (like TAIR).
# 3   feature The feature type name, like "gene" or "exon". In a well structured GFF file, all the children features always follow their parents in a single block (so all exons of a transcript are put after their parent "transcript" feature line and before any other parent transcript line). In GFF3, all features and their relationships should be compatible with the standards released by the Sequence Ontology Project.
# 4   start   Genomic start of the feature, with a 1-base offset. This is in contrast with other 0-offset half-open sequence formats, like BED.
# 5   end Genomic end of the feature, with a 1-base offset. This is the same end coordinate as it is in 0-offset half-open sequence formats, like BED.[citation needed]
# 6   score   Numeric value that generally indicates the confidence of the source in the annotated feature. A value of "." (a dot) is used to define a null value.
# 7   strand  Single character that indicates the strand of the feature; it can assume the values of "+" (positive, or 5'->3'), "-", (negative, or 3'->5'), "." (undetermined).
# 8   phase   phase of CDS features; it can be either one of 0, 1, 2 (for CDS features) or "." (for everything else). See the section below for a detailed explanation.
# 9   attributes  All the other information pertaining to this feature. The format, structure and content of this field is the one which varies the most between the three competing file formats.