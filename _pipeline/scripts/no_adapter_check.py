"""Stand-in for the porechop adapter check when it is off (adapter_check: false, and contigs).

Reads without the sequencing adapter (adapter-trimmed basecalls, public SRA/ENA reads) and assembled
contigs cannot show an adapter after the telomere, so porechop is skipped. This writes the two files
the telomere trimming step reads, for every anchored read:
  <output_dir>/<base>_adapter_trimming_check.tsv   same columns as compare_adapter_callers_dorado.py
                                                   writes; the adapter and tag columns are left blank
                                                   (not checked, which is not the same as False)
  <output_dir>/porechop/<base>_trimmed.fasta       the anchored reads, untrimmed

Usage: no_adapter_check.py <base_name> <anchor_set> <output_dir> <anchored_reads_fasta>
"""
import os
import shutil
import sys

import pandas as pd

COLUMNS = ['read_id', 'Repeat_Type', 'Adapter_After_Telomere', 'Tag_After_Telomere',
           'both_adapter_and_tag', 'either_adapter_or_tag', 'anchor_name',
           'dorado_telomere_adapter', 'dorado_centromere_adapter']
NOT_CHECKED = ['Adapter_After_Telomere', 'Tag_After_Telomere', 'both_adapter_and_tag', 'either_adapter_or_tag']


def adapter_table(top_matches):
    """The adapter-check table for every anchored read, with the adapter columns left blank."""
    df = pd.DataFrame({'read_id': top_matches['read_id'].astype(str),
                       'Repeat_Type': top_matches['repeat_type'],
                       'anchor_name': top_matches['anchor_name']})
    for col in NOT_CHECKED:
        df[col] = ''
    df['dorado_telomere_adapter'] = False
    df['dorado_centromere_adapter'] = False
    return df[COLUMNS]


def main():
    base_name, anchor_set, output_dir, anchored_fasta = sys.argv[1:5]
    print("Starting no_adapter_check.py (adapter check off)")
    top = pd.read_csv(os.path.join(output_dir, 'blast', f'top_matches_{base_name}_blasted_{anchor_set}.tsv'),
                      sep='\t')
    table = adapter_table(top)
    table.to_csv(os.path.join(output_dir, f'{base_name}_adapter_trimming_check.tsv'), sep='\t')
    porechop_dir = os.path.join(output_dir, 'porechop')
    os.makedirs(porechop_dir, exist_ok=True)
    shutil.copyfile(anchored_fasta, os.path.join(porechop_dir, f'{base_name}_trimmed.fasta'))
    print(f"  {len(table)} anchored reads passed on without an adapter check")


if __name__ == '__main__':
    main()
