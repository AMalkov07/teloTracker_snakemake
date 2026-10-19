#!/bin/bash
#$ -q UI,TELOMERE2
#$ -pe smp 48
#$ -j y
#$ -cwd
#$ -N yprime_pop
# Population-wide Y' count across ALL reads (not just anchored) for the 4 AM7575
# datasets. BLASTs the Y' probe against every read and tallies full-length
# (>=80bp) probe hits per read = Y' count. Captures Y'-bearing reads that never
# anchor (extrachromosomal Y' circles / concatemers) -- the strongest test for a
# Type I survivor's amplification.
#
# The 4 datasets run IN PARALLEL (the real speedup): the bottleneck is the
# single-threaded makeblastdb + I/O per dataset, not the near-instant blastn, so
# concurrency across datasets helps where extra blastn threads would not.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate consensus
export LANG=C LC_ALL=C

PROBE="_pipeline/references/y_prime_probe.fasta"
OUT="population_yprime_counts"
mkdir -p "$OUT"

process() {
    local tag="$1" base="$2"
    local READS="results/${base}/_pipeline/${base}.fasta"
    if [ ! -s "$READS" ]; then
        READS="${OUT}/${base}.fasta"
        samtools fasta -@ 8 "samples_dorado_basecalled/${base}.bam" > "$READS"
    fi
    if [ -s "${READS}.fai" ]; then wc -l < "${READS}.fai" > "${OUT}/${tag}_total.txt"
    else grep -c '^>' "$READS" > "${OUT}/${tag}_total.txt"; fi
    # Reads are the QUERY, probe is the tiny subject -> no max_target_seqs cap
    # (the probe-as-query direction silently capped at 500 reads). Each read is
    # searched against the 100bp probe; count >=80bp hits per read = Y' count.
    makeblastdb -in "$PROBE" -dbtype nucl -out "${OUT}/${tag}_probedb" >/dev/null 2>&1
    blastn -query "$READS" -db "${OUT}/${tag}_probedb" -task blastn \
        -outfmt "6 qseqid length" -num_threads 10 2>/dev/null \
        | awk '$2>=80{c[$1]++} END{for(r in c) print c[r]}' > "${OUT}/${tag}_percounts.txt"
    rm -f "${OUT}/${tag}_probedb".*
    echo "[$tag] done: $(cat ${OUT}/${tag}_total.txt) total reads, $(wc -l < ${OUT}/${tag}_percounts.txt) with >=1 Y'"
}

# launch all 4 concurrently
process A3A_D0 dorado_AM7575_A3A-3-D0_072326 &
process EV_D0  dorado_AM7575_EV-3-D0_072226 &
process A3A_4c dorado_AM7575_A3A-3-4c &
process EV_4c  dorado_AM7575_EV-3-4c &
wait
echo "all 4 datasets processed"

python3 - <<'PY'
import collections
OUT="population_yprime_counts"
tags=["A3A_D0","EV_D0","A3A_4c","EV_4c"]
hist={}; maxc=0
for t in tags:
    total=int(open(f"{OUT}/{t}_total.txt").read().strip())
    counts=[int(x) for x in open(f"{OUT}/{t}_percounts.txt")]
    d=collections.Counter(counts)
    d[0]=total-sum(d.values())
    hist[t]=d; maxc=max(maxc, max(d) if d else 0)
with open(f"{OUT}/yprime_population_distribution.tsv","w") as f:
    f.write("yprime_count\t"+"\t".join(tags)+"\n")
    for k in range(0,maxc+1):
        f.write(f"{k}\t"+"\t".join(str(hist[t].get(k,0)) for t in tags)+"\n")
print("\nwrote population_yprime_counts/yprime_population_distribution.tsv\n")
print("=== % of ALL reads (population-wide) with >=N Y' ===")
print("N    "+"  ".join(f"{t:>8}" for t in tags))
for N in [1,3,6,10,15,20]:
    row=[]
    for t in tags:
        tot=sum(hist[t].values()); c=sum(v for k,v in hist[t].items() if k>=N)
        row.append(f"{100*c/tot:>7.3f}%")
    print(f">={N:<3}"+"  ".join(row))
PY
echo "DONE"
