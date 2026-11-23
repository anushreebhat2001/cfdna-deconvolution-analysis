"""Number every CpG in hg38 the way Loyfer/wgbstools does, and check against the authors' file.
Output: BED (chr, start, end, cpg_index) with start = 0-based position of the C."""
import sys, gzip, re, collections

fasta, netanel, out = sys.argv[1], sys.argv[2], sys.argv[3]

# expected chromosome order, count and first index, from the authors' CpG-number file
cnt, first = collections.OrderedDict(), {}
with gzip.open(netanel, "rt") as f:
    for line in f:
        c, i = line.split()
        if c not in cnt:
            cnt[c] = 0
            first[c] = int(i)
        cnt[c] += 1
print("authors' file:", len(cnt), "chromosomes,", sum(cnt.values()), "CpGs")

# check the numbering is continuous across chromosomes in that order
run = 1
for c in cnt:
    if first[c] != run:
        print("WARNING: index jump at", c, "expected", run, "got", first[c])
    run += cnt[c]

ours = {}
fo = open(out, "w")
def process(name, parts):
    if name not in cnt:
        return
    seq = "".join(parts).upper()
    pos = [m.start() for m in re.finditer("CG", seq)]
    ours[name] = len(pos)
    status = "OK" if len(pos) == cnt[name] else "MISMATCH"
    print(f"{name}: ours {len(pos)}  authors {cnt[name]}  {status}", flush=True)
    if status == "OK":
        base = first[name]
        for k, p in enumerate(pos):
            fo.write(f"{name}\t{p}\t{p+1}\t{base+k}\n")

name, parts = None, []
with gzip.open(fasta, "rt") as f:
    for line in f:
        if line.startswith(">"):
            if name:
                process(name, parts)
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.rstrip())
    if name:
        process(name, parts)
fo.close()

missing = [c for c in cnt if c not in ours]
bad = [c for c in ours if ours[c] != cnt[c]]
print("chromosomes missing from fasta:", missing)
print("chromosomes with count mismatch:", bad)
print("ALL MATCH" if not missing and not bad else "NOT ALL MATCH")
