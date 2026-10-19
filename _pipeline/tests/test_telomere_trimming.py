"""Unit tests for the telomere tract search (fine_telomere_trimming.trim_for_telomere_repeat).
Run: python _pipeline/tests/test_telomere_trimming.py

A read with no further repeat after its tract, or none at all, used to raise inside the search and
was silently dropped from every output. It must now come back as a tract or as "no telomere" (None).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import fine_telomere_trimming as ftt

AC_TELO = 'CCCACACCACACCCACACCCACACCACACCCA' * 8          # 256 bp of C1-3A repeat
AC_FLANK = 'GATTAGTTGTGAGTTGATTGTGATGAGTTGTAGTGAGTTG' * 10  # no C: cannot extend an AC tract
TG_FLANK = 'CAATCACTTCAACTCATACTCAATTCACTAACTCATCAAC' * 10  # no G: cannot extend a TG tract


def revcomp(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def test_ac_read_tract_and_offset():
    seq, length, _, offset = ftt.trim_for_telomere_repeat('TTAG' + AC_TELO + AC_FLANK, 'AC')
    assert (length, offset) == (len(AC_TELO), 4) and seq == AC_TELO


def test_tg_read_tract_and_offset():
    read = TG_FLANK + revcomp(AC_TELO) + 'CTAA'
    seq, length, _, offset = ftt.trim_for_telomere_repeat(read, 'TG')
    assert (length, offset) == (len(AC_TELO), 4) and seq == revcomp(AC_TELO)


def test_all_telomere_read_is_kept():
    # nothing follows the tract: it runs to the read end instead of the read being dropped
    _, length, _, offset = ftt.trim_for_telomere_repeat(AC_TELO, 'AC')
    assert (length, offset) == (len(AC_TELO), 0)
    _, length, _, _ = ftt.trim_for_telomere_repeat(revcomp(AC_TELO), 'TG')
    assert length == len(AC_TELO)


def test_read_without_repeat_is_no_telomere():
    assert ftt.trim_for_telomere_repeat('G' * 500, 'AC') is None
    assert ftt.trim_for_telomere_repeat('', 'AC') is None
    assert ftt.trim_for_telomere_repeat('', 'TG') is None


def test_seed_too_far_from_read_end_is_no_telomere():
    assert ftt.trim_for_telomere_repeat('G' * 100 + AC_TELO, 'AC') is None


def test_lowercase_matches_uppercase():
    read = 'TTAG' + AC_TELO + AC_FLANK
    assert ftt.trim_for_telomere_repeat(read.lower(), 'AC')[1] == ftt.trim_for_telomere_repeat(read, 'AC')[1]


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
