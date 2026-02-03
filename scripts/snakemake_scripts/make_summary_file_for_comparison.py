import pandas as pd
from Bio import SeqIO
import sys

def fasta_to_length_tsv(fasta_path, output_tsv):
    records = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        records.append({
            "read_id": record.id,
            "sequence_length_template": len(record.seq)
        })

    df = pd.DataFrame(records)
    df.to_csv(output_tsv, sep="\t", index=False)



###################

#base_name = sys.argv[1]

#fasta_file = f'{base_name}/{base_name}.fasta'
fasta_file = sys.argv[1]
#summary_tsv = f'{base_name}/{base_name}_sequencing_summary.tsv'
summary_tsv = sys.argv[2]

fasta_to_length_tsv(fasta_file, summary_tsv)
