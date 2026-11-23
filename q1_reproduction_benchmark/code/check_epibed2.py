import sys, re, bisect, collections, pysam

epibed, atlas, regions = sys.argv[1], sys.argv[2], sys.argv[3]

cpgs = set()
with open(atlas) as f:
    next(f)
    for line in f:
        c, s = line.split("\t", 2)[:2]
        cpgs.add((c, int(s)))

regs = []
for l in open(regions):
    if l.strip():
        c, s, e = l.split()[:3]
        regs.append((c, int(s), int(e)))
byc = collections.defaultdict(list)
for c, s, e in regs:
    byc[c].append((s, e))
starts = {}
for c in byc:
    byc[c].sort()
    starts[c] = [s for s, e in byc[c]]

def in_region(c, p):
    v = byc.get(c)
    if not v:
        return False
    i = bisect.bisect_right(starts[c], p) - 1
    return i >= 0 and v[i][0] <= p < v[i][1]

tb = pysam.TabixFile(epibed)
TOK = re.compile(r"([A-Za-z])(\d*)")
tot, hit = collections.Counter(), collections.Counter()
seen, examples = set(), []

for c, s, e in regs:
    try:
        rows = tb.fetch(c, s, e)
    except ValueError:
        continue
    for line in rows:
        f = line.rstrip("\n").split("\t")
        if f[3] in seen:
            continue
        seen.add(f[3])
        pos, inside = int(f[1]), []
        for ch, num in TOK.findall(f[6]):
            n = int(num) if num else 1
            if ch == "i":
                continue
            for k in range(n):
                p = pos + k
                if in_region(c, p):
                    tot[ch] += 1
                    inside.append((p, ch))
                    for sh in (-1, 0, 1):
                        if (c, p + sh) in cpgs:
                            hit[(ch, sh)] += 1
            pos += n
        if len(examples) < 4 and any(ch != "x" for p, ch in inside):
            examples.append((c, s, e, f[1], f[2], f[5], f[6], [(p, ch) for p, ch in inside if ch != "x"]))

print("reads overlapping markers:", len(seen))
print("positions INSIDE marker regions, by code (share landing on an atlas CpG row at shift -1 / 0 / +1):")
for ch in sorted(tot):
    print(f"  {ch}: {tot[ch]:>7} positions | " +
          " / ".join(f"{100*hit[(ch, sh)]/tot[ch]:5.1f}%" for sh in (-1, 0, 1)))

for c, s, e, rs, re_, strand, rle, calls in examples:
    region_cpgs = sorted(p for (cc, p) in cpgs if cc == c and s <= p < e)
    print(f"\nregion {c}:{s}-{e}  atlas CpG rows: {region_cpgs}")
    print(f"  read {rs}-{re_} strand {strand}")
    print(f"  non-x codes inside region: {calls}")
