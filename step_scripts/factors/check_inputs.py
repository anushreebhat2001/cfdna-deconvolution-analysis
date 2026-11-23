import contextlib
import gzip
import io
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "step0"))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import cut_strings
from score_pair import ATLAS
sys.path.insert(0, HERE)
from make_inputs import BREAST, CONDITIONS, OUT

REGIONS = "data/ref/31CT_U250_Netanel_merged_regions_file.bed"
TMP = tempfile.mkdtemp()
ok = True


def report(name, passed, detail=""):
    global ok
    ok &= bool(passed)
    print(f"{name}: {'pass' if passed else 'FAIL'} {detail}")


def cut(src, *options):
    out = os.path.join(TMP, "out.pat")
    sys.argv = ["cut_strings.py", src, REGIONS, out, *options]
    with contextlib.redirect_stdout(io.StringIO()):
        cut_strings.main()
    return open(out).read()


def beads(text):
    out = []
    for line in text.splitlines():
        _, s, p, n = line.split("\t")
        out += [(int(s) + i, c, n) for i, c in enumerate(p) if c != "."]
    return sorted(out)


print("A  modified atlases")
orig = [line.rstrip("\n").split("\t") for line in open(ATLAS["u250"])]
for tag, value in (("off", "5"), ("on", "95")):
    new = [line.rstrip("\n").split("\t") for line in open(f"{OUT}/atlas/u250_plus_{tag}.txt")]
    report(f"A u250_plus_{tag}: 31 original types unchanged at {len(orig) - 1} CpGs",
           len(new) == len(orig) and all(n[:-2] == o for n, o in zip(new, orig)))
    report(f"A u250_plus_{tag}: added column is {value}/100 at every CpG",
           all(n[-2:] == [value, "100"] for n in new[1:]), f"(header {new[0][-2:]})")
orig = [line.rstrip("\n").split("\t") for line in open(ATLAS["max"])]
new = [line.rstrip("\n").split("\t") for line in open(f"{OUT}/atlas/max_minus_megak.txt")]
report("A max_minus_megak: the 7 other types unchanged, megakaryocyte column gone",
       len(new) == len(orig) and all(n == o[:-2] for n, o in zip(new, orig)) and "Megakaryocytes_METH" not in new[0],
       f"(types left: {[h[:-5] for h in new[0] if h.endswith('_METH')]})")

print("\nB  the 46 files of T7")
meta = [line.rstrip("\n").split("\t") for line in list(open("meta/samples.tsv"))[1:]]
files = [(f[0], "results/step2") for f in meta] + [(s, "results/breast") for s in BREAST]
same_whole = same_cut = same_beads = 0
for s, d in files:
    whole, single = cut(f"data/marker_pat/{s}.pat"), cut(f"data/marker_pat/{s}.pat", "--size", "1")
    same_whole += whole == gzip.open(f"{d}/{s}/whole.pat.gz", "rt").read()
    same_cut += single == gzip.open(f"{d}/{s}/cut.pat.gz", "rt").read()
    same_beads += beads(whole) == beads(single)
report("B whole strings cut again equal the fitted ones", same_whole == len(files), f"({same_whole} of {len(files)})")
report("B single beads cut again equal the fitted ones", same_cut == len(files), f"({same_cut} of {len(files)})")
report("B whole and single hold the same beads", same_beads == len(files), f"({same_beads} of {len(files)})")

print("\nC  breast dilutions")
pool_main = [f[0] for f in meta if f[1] == "control" and f[6] == "yes"]
pool_breast = [f[0] for f in (line.rstrip("\n").split("\t") for line in list(open(f"{OUT}/meta_breast.tsv"))[1:])
               if f[1] == "control" and f[6] == "yes"]
report("C same 20 pooled controls as the cancer dilutions", pool_main == pool_breast and len(pool_main) == 20)


def shares(path):
    out = {}
    for line in list(open(path))[1:]:
        s, pct, n, mc, mh, kc, kh = line.rstrip("\n").split("\t")
        out.setdefault(pct, []).append((int(mc) / (int(mc) + int(mh)), int(kc) / (int(kc) + int(kh)), int(kc) + int(kh)))
    return out


print("level\tbreast: read share from sample\tbreast: string share from sample (range)\tbreast: strings per file (range)\t"
      "cancer dilutions: string share (range)\tcancer dilutions: strings per file (range)")
b, c = shares("data/mix_breast/mix_log.tsv"), shares("data/mix/mix_log.tsv")
for pct in b:
    x, y = np.array(b[pct]), np.array(c[pct])
    print(f"{pct}%\t{x[:, 0].mean():.4f}\t{x[:, 1].min():.3f}-{x[:, 1].max():.3f}\t{int(x[:, 2].min())}-{int(x[:, 2].max())}\t"
          f"{y[:, 1].min():.3f}-{y[:, 1].max():.3f}\t{int(y[:, 2].min())}-{int(y[:, 2].max())}")

print("\nD  simulation")
controls = [f[0] for f in meta if f[1] == "control"]
same = sum(open(f"{OUT}/sim/none/{s}.pat").read() == open(f"results/d1/sim/{s}.pat").read() for s in controls)
report("D condition 'none' reproduces the stored D1 files", same == len(controls), f"({same} of {len(controls)})")
print("condition\tstrings\tbeads\tstrings whose values changed\tbeads on->off\tbeads off->on\tall-off strings among 2+ bead strings")
for cond in CONDITIONS:
    shape = n_str = n_bead = ch_str = on_off = off_on = alloff = multi = 0
    for s in controls:
        base = [line.rstrip("\n").split("\t") for line in open(f"{OUT}/sim/none/{s}.pat")]
        new = [line.rstrip("\n").split("\t") for line in open(f"{OUT}/sim/{cond}/{s}.pat")]
        shape += len(base) == len(new) and all(
            x[:2] == y[:2] and x[3] == y[3] and [c == "." for c in x[2]] == [c == "." for c in y[2]] for x, y in zip(base, new))
        for x, y in zip(base, new):
            n_str += 1
            n_bead += len(x[2]) - x[2].count(".")
            ch_str += x[2] != y[2]
            on_off += sum(a == "C" and b == "T" for a, b in zip(x[2], y[2]))
            off_on += sum(a == "T" and b == "C" for a, b in zip(x[2], y[2]))
            v = y[2].replace(".", "")
            if len(v) >= 2:
                multi += 1
                alloff += set(v) == {"T"}
    report(f"D {cond}: places, counts and bead positions unchanged", shape == len(controls))
    print(f"{cond}\t{n_str}\t{n_bead}\t{ch_str} ({100 * ch_str / n_str:.2f}%)\t{on_off} ({100 * on_off / n_bead:.2f}%)\t"
          f"{off_on} ({100 * off_on / n_bead:.2f}%)\t{100 * alloff / multi:.2f}%")
print("\nALL CHECKS PASSED" if ok else "\nSOME CHECK FAILED")
