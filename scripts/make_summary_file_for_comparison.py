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

base_name = sys.argv[1]

fasta_file = f'results/{base_name}/{base_name}.fasta'
summary_tsv = f'results/{base_name}/{base_name}_sequencing_summary.tsv'

fasta_to_length_tsv(fasta_file, summary_tsv)
