import collections
import contextlib
import io
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "step0"))
import cut_strings

REGIONS = "data/ref/31CT_U250_Netanel_merged_regions_file.bed"
TMP = tempfile.mkdtemp()


def cut(src, *options):
    out = os.path.join(TMP, "out.pat")
    sys.argv = ["cut_strings.py", src, REGIONS, out, *options]
    with contextlib.redirect_stdout(io.StringIO()):
        cut_strings.main()
    return [(int(f[1]), f[2]) for f in (line.split("\t") for line in open(out))]


def beads(lines):
    return sorted((s + i, c) for s, p in lines for i, c in enumerate(p) if c != ".")


def n_beads(pattern):
    return len(pattern) - pattern.count(".")


def main():
    inputs = subprocess.run(["bash", "step_scripts/common/inputs.sh"], capture_output=True, text=True).stdout.split()
    fail = collections.Counter()
    targets, short, smallest = [], collections.defaultdict(list), collections.Counter()
    for src in inputs:
        everything = cut(src)
        full = {b: cut(src, "--bucket", b) for b in ("1", "2", "3", "4+")}
        for b, lines in full.items():
            ok = all(n_beads(p) >= 4 if b == "4+" else n_beads(p) == int(b) for _, p in lines)
            fail["U1"] += not ok
        for b in ("2", "3"):
            single = cut(src, "--bucket", b, "--size", "1")
            fail["U2"] += not (beads(single) == beads(full[b]) and all(len(p) == 1 for _, p in single))
        fail["U3"] += beads(everything) != sorted(sum((beads(v) for v in full.values()), []))

        size = {b: len(beads(full[b])) for b in ("2", "3", "4+")}
        target = min(size.values())
        targets.append(target)
        smallest[min(size, key=size.get)] += 1
        for b in ("2", "3", "4+"):
            draw1 = cut(src, "--bucket", b, "--match", "1")
            draw2 = cut(src, "--bucket", b, "--match", "2")
            again = cut(src, "--bucket", b, "--match", "1")
            single = cut(src, "--bucket", b, "--match", "1", "--size", "1")
            fail["M1"] += bool(collections.Counter(draw1) - collections.Counter(full[b]))
            got = len(beads(draw1))
            longest = max(n_beads(p) for _, p in full[b])
            fail["M2"] += not (target - longest < got <= target)
            short[b].append(target - got)
            fail["M3"] += not (beads(single) == beads(draw1) and all(len(p) == 1 for _, p in single))
            thinned = size[b] > target
            same12 = sorted(draw1) == sorted(draw2)
            fail["M4"] += again != draw1 or (thinned and same12) or (not thinned and not same12)

    print(f"inputs checked: {len(inputs)}")
    for t in ("U1", "U2", "U3", "M1", "M2", "M3", "M4"):
        print(f"{t}: {'pass' if fail[t] == 0 else f'FAIL on {fail[t]} cases'}")
    t = np.array(targets)
    print(f"matched bead total per input: smallest {t.min()}, median {int(np.median(t))}, largest {t.max()}")
    print("which bucket is the smallest (sets the total):", dict(smallest))
    for b in ("2", "3", "4+"):
        print(f"bucket {b}: beads short of the total after thinning: at most {max(short[b])}")


if __name__ == "__main__":
    main()
