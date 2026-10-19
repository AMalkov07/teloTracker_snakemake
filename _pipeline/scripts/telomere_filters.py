"""Which anchored reads are 'qualifying': reads whose telomere was found and measured.

One definition for every script that draws or counts them (y_prime_analysis.py, single_sample_plots.py,
read_summary.py, make_sample_report.py): a telomere repeat of at least MIN_TELOMERE_REPEAT bp and, when
the adapter check ran, the sequencing adapter right after it, which shows the read reached the end of
the molecule. Without the adapter check (adapter-trimmed or public reads, contigs) only the repeat
length counts; the telomere trimming step already requires the repeat to reach the read end.
"""
import pandas as pd

MIN_TELOMERE_REPEAT = 30


def adapter_found(values):
    """True where the Adapter_After_Telomere column says True (bool, or 'True' / '1' as text)."""
    return values.astype(str).str.strip().str.lower().isin(['true', '1'])


def qualifying_mask(df, require_adapter=True):
    """Boolean mask over a per-read table with repeat_length and Adapter_After_Telomere columns."""
    mask = pd.to_numeric(df['repeat_length'], errors='coerce') >= MIN_TELOMERE_REPEAT
    if require_adapter:
        mask &= adapter_found(df['Adapter_After_Telomere'])
    return mask
