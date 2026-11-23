import gzip
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stage1_checks as s1
import stage2_readout as r


def raw_calls():
    print("1. raw calls per group: median over samples (min-max)")
    print("group\tsamples\tshare of beads called methylated\tall-on share of 2+ bead strings\tall-off share of 2+ bead strings")
    for g, sel in r.GROUPS.items():
        fr, on, off = [], [], []
        for i in sel:
            d = "results/breast" if i in r.BREAST else "results/step2"
            C = T = n2 = allon = alloff = 0
            with gzip.open(f"{d}/{i}/whole.pat.gz", "rt") as f:
                for line in f:
                    p = line.split("\t")
                    k, c, t = int(p[3]), p[2].count("C"), p[2].count("T")
                    C += c * k
                    T += t * k
                    if c + t >= 2:
                        n2 += k
                        allon += k * (t == 0)
                        alloff += k * (c == 0)
            fr.append(C / (C + T))
            on.append(allon / n2)
            off.append(alloff / n2)
        f3 = lambda v: f"{np.median(v):.3f} ({min(v):.3f}-{max(v):.3f})"
        print(f"{g}\t{len(sel)}\t{f3(fr)}\t{f3(on)}\t{f3(off)}")


def t7_off_by_sample():
    print("\n2. T7 with the off column: epithelial score (16 types) per sample, and how many of the 22 controls it is above")
    epi = np.concatenate([s1.MASK["epithelial16"], [False]])
    ids = [m["id"] for m in r.META] + r.BREAST
    tf = {m["id"]: m["tf"] for m in r.META}
    pat = {m["id"]: m["patient"] for m in r.META}
    folder = lambda i: f"{r.O}/t7_off"
    for clip in ("unclipped", "clipped"):
        sc = r.scores(folder, epi, ids, clip)
        off = r.scores(folder, np.arange(32) == 31, ids, clip)
        ctrl = np.array([sc[m["id"]] for m in r.CTRL])
        print(f"\n{clip}; controls min/median/max: whole {ctrl[:, 0].min():.3f}/{np.median(ctrl[:, 0]):.3f}/{ctrl[:, 0].max():.3f}, "
              f"single {ctrl[:, 1].min():.3f}/{np.median(ctrl[:, 1]):.3f}/{ctrl[:, 1].max():.3f}")
        print("sample\tpatient\ttumor share\twhole\tsingle\tcontrols below, whole\tcontrols below, single\toff column whole\toff column single")
        rows = r.BREAST + sorted((m["id"] for m in r.CANC if m["tf"] < r.LOW), key=tf.get)
        for i in rows:
            w, c = sc[i]
            t = f"{tf[i]:.3f}" if i in tf else "0-0.03"
            print(f"{i}\t{pat.get(i, 'breast')}\t{t}\t{w:.3f}\t{c:.3f}\t{(w > ctrl[:, 0]).sum()}\t{(c > ctrl[:, 1]).sum()}\t{off[i][0]:.3f}\t{off[i][1]:.3f}")


if __name__ == "__main__":
    print("EXPLORATORY follow-up to stage 2, not in PLAN.md (made after seeing results/factors/stage2.txt)\n")
    raw_calls()
    t7_off_by_sample()
