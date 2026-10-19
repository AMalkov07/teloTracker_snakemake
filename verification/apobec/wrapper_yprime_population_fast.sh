#!/bin/bash
#$ -q UI,TELOMERE2
#$ -pe smp 24
#$ -j y
#$ -cwd
#$ -N yprime_pop_fast
# FASTER population-wide Y' count, run IN PARALLEL with the reads-as-query job.
# Direction: probe is the QUERY (1 query), reads are the DB. This is far fewer
# blastn queries than the reads-as-query direction, but the probe-as-query
# direction previously CAPPED at 500 reads via the default -max_target_seqs.
# Fix = raise the cap to 5,000,000 so every read with a probe hit is reported.
# Writes to a SEPARATE output dir so it never touches the running job's files.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate consensus
export LANG=C LC_ALL=C

PROBE="_pipeline/references/y_prime_probe.fasta"
OUT="population_yprime_fast"
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
    # reads as DB (built once); probe as the single query. High max_target_seqs
    # removes the 500-read cap. Count >=80bp HSPs per read (sseqid) = Y' count.
    makeblastdb -in "$READS" -dbtype nucl -out "${OUT}/${tag}_readsdb" >/dev/null 2>&1
    blastn -query "$PROBE" -db "${OUT}/${tag}_readsdb" -task blastn \
        -max_target_seqs 5000000 -outfmt "6 sseqid length" -num_threads 6 2>/dev/null \
        | awk '$2>=80{c[$1]++} END{for(r in c) print c[r]}' > "${OUT}/${tag}_percounts.txt"
    rm -f "${OUT}/${tag}_readsdb".*
    echo "[$tag] done: $(cat ${OUT}/${tag}_total.txt) total reads, $(wc -l < ${OUT}/${tag}_percounts.txt) with >=1 Y'"
}

process A3A_D0 dorado_AM7575_A3A-3-D0_072326 &
process EV_D0  dorado_AM7575_EV-3-D0_072226 &
process A3A_4c dorado_AM7575_A3A-3-4c &
process EV_4c  dorado_AM7575_EV-3-4c &
wait
echo "all 4 datasets processed"

python3 - <<'PY'
import collections
OUT="population_yprime_fast"
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
print("\nwrote population_yprime_fast/yprime_population_distribution.tsv\n")
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
