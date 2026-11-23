import sys, re, collections, pysam

epibed, atlas, regions = sys.argv[1], sys.argv[2], sys.argv[3]

# atlas CpG positions (0-based C position, as in the atlas START column)
cpgs = set()
with open(atlas) as f:
    next(f)
    for line in f:
        c, s = line.split("\t", 2)[:2]
        cpgs.add((c, int(s)))
print("atlas CpGs:", len(cpgs))

regs = [l.split()[:3] for l in open(regions) if l.strip()]
print("marker regions:", len(regs))

tb = pysam.TabixFile(epibed)
TOK = re.compile(r"([A-Za-z])(\d*)")
chars = collections.Counter()
calls = collections.Counter()           # M/U calls per strand
hits = collections.Counter()            # (strand, shift) -> calls landing on an atlas CpG
cpg_per_read = collections.Counter()    # number of M/U calls inside marker CpGs, per read
seen, nreads, len_mismatch = set(), 0, 0

for c, s, e in regs:
    try:
        rows = tb.fetch(c, int(s), int(e))
    except ValueError:
        continue
    for line in rows:
        f = line.rstrip("\n").split("\t")
        if f[3] in seen:
            continue
        seen.add(f[3]); nreads += 1
        start, end, strand, rle = int(f[1]), int(f[2]), f[5], f[6]
        pos, in_marker = start, 0
        for ch, num in TOK.findall(rle):
            n = int(num) if num else 1
            chars[ch] += n
            if ch == "i":              # insertion: uses no reference positions
                continue
            if ch in "MU":
                for k in range(n):
                    p = pos + k
                    calls[strand] += 1
                    for sh in (-1, 0, 1):
                        if (f[0], p + sh) in cpgs:
                            hits[(strand, sh)] += 1
                    if (f[0], p) in cpgs or (f[0], p - 1) in cpgs:
                        in_marker += 1
            pos += n
        if pos != end:
            len_mismatch += 1
        cpg_per_read[in_marker] += 1

print("reads overlapping markers:", nreads)
print("reads whose decoded length != end-start:", len_mismatch)
print("RLE character totals:", dict(chars))
for strand in sorted(calls):
    print(f"strand {strand}: {calls[strand]} M/U calls")
    for sh in (-1, 0, 1):
        print(f"   shift {sh:+d}: {hits[(strand, sh)]} hits ({100*hits[(strand, sh)]/max(calls[strand],1):.1f}%)")
print("marker CpGs per read (count: n_reads):", dict(sorted(cpg_per_read.items())))
