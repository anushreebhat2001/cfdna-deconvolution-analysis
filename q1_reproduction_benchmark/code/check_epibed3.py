import sys, re, gzip, collections, pysam

cpg_index = sys.argv[1]
files = sys.argv[2:]
N = 200000                        # first 200k reads of each file (all on chr1)
TOK = re.compile(r"([A-Za-z])(\d*)")

cpgC = set()
for row in pysam.TabixFile(cpg_index).fetch("chr1"):
    cpgC.add(int(row.split("\t")[1]))
print("chr1 CpGs loaded:", len(cpgC))

for path in files:
    tot, onC, onG = collections.Counter(), collections.Counter(), collections.Counter()
    calls_per_read = collections.Counter()
    indel, ref_len = 0, 0
    with gzip.open(path, "rt") as f:
        for n, line in enumerate(f):
            if n >= N:
                break
            fld = line.rstrip("\n").split("\t")
            if fld[0] != "chr1":
                break
            pos, ncall = int(fld[1]), 0
            for ch, num in TOK.findall(fld[6]):
                k = int(num) if num else 1
                if ch == "i":
                    indel += k
                    continue
                if ch == "d":
                    indel += k
                ref_len += k
                if ch != "x":
                    for j in range(k):
                        p = pos + j
                        tot[ch] += 1
                        if p in cpgC: onC[ch] += 1
                        elif p - 1 in cpgC: onG[ch] += 1
                if ch in "MU":
                    ncall += k
                pos += k
            calls_per_read[min(ncall, 6)] += 1
    print(f"\n=== {path.split('/')[-1]}  ({n} reads)")
    print(f"  indel rate: {100*indel/ref_len:.2f}% of bases")
    for ch in sorted(tot):
        print(f"  code {ch}: {tot[ch]:>8} | on CpG C {100*onC[ch]/tot[ch]:5.1f}% | on CpG G {100*onG[ch]/tot[ch]:5.1f}%")
    s = sum(calls_per_read.values())
    print("  M/U calls per read (6 = 6+):", {k: f"{100*v/s:.1f}%" for k, v in sorted(calls_per_read.items())})
