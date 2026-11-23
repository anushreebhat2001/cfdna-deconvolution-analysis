import multiprocessing
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from fit_linked import ISH_ITERATIONS, ISH_STOP, LINKED_ITERATIONS, ish_start, load
from linked_model import LinkedISH
from score_pair import ATLAS, EPI

REGIONS = "data/ref/31CT_U250_Netanel_merged_regions_file.bed"
NAMES = [h[:-5] for h in open(ATLAS["u250"]).readline().rstrip("\n").split("\t") if h.endswith("_METH")]
EPI_MASK = np.array([n in EPI["u250"] for n in NAMES])
TUMOR = NAMES.index("Colon-Ep")
META = [line.rstrip("\n").split("\t") for line in list(open("meta/samples.tsv"))[1:]]
RHOS, SHARES, REPS = (0.0, 0.1, 0.3, 0.6), (0.0, 0.02, 0.10), 20
SHARES_STACKED, STACKS = (0.0, 0.10), 8


def six_files():
    ctrl = {b: next(f[0] for f in META if f[1] == "control" and f[3] == b) for b in ("HU012", "HU083")}
    crc = sorted((f for f in META if f[1] == "crc"), key=lambda f: float(f[5]))
    return [("control HU012", f"results/step2/{ctrl['HU012']}"), ("control HU083", f"results/step2/{ctrl['HU083']}"),
            ("cancer, highest tumor share", f"results/step2/{crc[-1][0]}"),
            ("cancer, lowest tumor share", f"results/step2/{crc[0][0]}"),
            ("3.125% dilution", f"results/step2/{crc[-1][0]}__d3.125"),
            ("D1 simulated control", f"results/d1/{ctrl['HU012']}")]


def stored(prefix, clip):
    return np.loadtxt(f"{prefix}.ish_{clip}.txt").ravel()


def real_file_checks(job):
    label, folder = job
    out = [f"\n{label}: {folder}"]
    eng = LinkedISH(*load(f"{folder}/whole.pat.gz", ATLAS["u250"], REGIONS))
    start = ish_start(eng.t)
    a0, _, it0, _, _ = eng.fit(start, 0.0, False, ISH_ITERATIONS, ISH_STOP)
    d = [np.abs(np.round(a0 * 100, 4) - stored(f"{folder}/whole", "unclipped")).max(),
         np.abs(np.round(eng.clipped(a0, 0.0) * 100, 4) - stored(f"{folder}/whole", "clipped")).max()]
    ok1 = max(d) <= 0.0001
    out.append(f"  1  rho = 0 against stored CelFiE-ISH: largest difference {d[0]:.4f} (unclipped), {d[1]:.4f} (clipped) "
               f"percent points, {it0} EM rounds  {'PASS' if ok1 else 'FAIL'}")
    cut = LinkedISH(*load(f"{folder}/cut.pat.gz", ATLAS["u250"], REGIONS))
    d2 = []
    for rho in (0.3, 0.9):
        a, _, _, _, _ = cut.fit(ish_start(cut.t), rho, False, ISH_ITERATIONS, ISH_STOP)
        d2 += [np.abs(np.round(a * 100, 4) - stored(f"{folder}/cut", "unclipped")).max(),
               np.abs(np.round(cut.clipped(a, rho) * 100, 4) - stored(f"{folder}/cut", "clipped")).max()]
    ok2 = max(d2) <= 0.0001 and not cut.multi.any()
    out.append(f"  2  single-bead file with rho fixed at 0.3 and 0.9: largest difference to stored {max(d2):.4f}, "
               f"strings with 2+ beads {int(cut.multi.sum())}  {'PASS' if ok2 else 'FAIL'}")
    alpha, rho, it, conv, trace = eng.fit(start, 0.5, True, LINKED_ITERATIONS, ISH_STOP, track=True)
    drop = min(np.diff(trace))
    ll0, ll1 = eng.loglik(a0, 0.0), eng.loglik(alpha, rho)
    ok3 = drop >= -1e-9 * abs(ll1) and ll1 >= ll0 - 1e-6 and conv
    out.append(f"  3  linked fit: rho {rho:.4f}, {it} EM rounds, converged {'yes' if conv else 'NO'}; largest drop of the "
               f"likelihood between rounds {min(drop, 0):.2e}; log-likelihood CelFiE-ISH {ll0:.3f}, linked {ll1:.3f} "
               f"({ll1 - ll0:+.3f})  {'PASS' if ok3 else 'FAIL'}")
    fits = []
    for seed in (1, 2):
        rng = np.random.default_rng(seed)
        s = rng.uniform(size=eng.t)
        for r0 in (0.02, 0.5, 0.98):
            a, r, _, c, _ = eng.fit(s / s.sum(), r0, True, LINKED_ITERATIONS, ISH_STOP)
            fits.append((a, r, c))
    spread_a = max(np.abs(a - alpha).max() for a, _, _ in fits) * 100
    spread_r = max(abs(r - rho) for _, r, _ in fits)
    ok4 = spread_a <= 0.01 and spread_r <= 0.0005 and all(c for _, _, c in fits)
    out.append(f"  4  six other starts: largest change of any cell-type share {spread_a:.4f} percent points, of rho "
               f"{spread_r:.5f}  {'PASS' if ok4 else 'FAIL'}")
    grid = np.round(np.arange(0, 0.96, 0.05), 2)
    prof = [eng.loglik(eng.fit(start, g, False, LINKED_ITERATIONS, ISH_STOP)[0], g) for g in grid]
    best = int(np.argmax(prof))
    step = np.sign(np.round(np.diff(prof), 6))
    one_peak = not (np.diff(step[step != 0]) > 0).any()
    ok5 = abs(grid[best] - rho) <= 0.05 and ll1 >= max(prof) - 1e-6 and one_peak
    out.append(f"  5  likelihood over fixed rho (0 to 0.95 in steps of 0.05): highest at {grid[best]:.2f}, fitted rho "
               f"{rho:.4f}, fitted likelihood minus best on the grid {ll1 - max(prof):+.4f}, one peak "
               f"{'yes' if one_peak else 'NO'}  {'PASS' if ok5 else 'FAIL'}")
    return "\n".join(out), all((ok1, ok2, ok3, ok4, ok5))


_SIM = {}


def sim_setup():
    ctrl = next(f[0] for f in META if f[1] == "control")
    eng = LinkedISH(*load(f"results/step2/{ctrl}/whole.pat.gz", ATLAS["u250"], REGIONS))
    healthy = np.mean([stored(f"results/step2/{f[0]}/whole", "unclipped") for f in META if f[1] == "control"], axis=0)
    healthy[EPI_MASK] = 0
    healthy /= healthy.sum()
    rates = []
    for b in eng.beta:
        known = ~np.isnan(b)
        fill = np.where(known, b, 0).sum(0) / np.maximum(known.sum(0), 1)
        rates.append(np.where(known, b, fill))
    _SIM.update(present=[x != 0 for x in eng.x], beta=eng.beta, rates=rates, healthy=healthy, ctrl=ctrl)


def make_file(rho, share, rng):
    mix = _SIM["healthy"] * (1 - share)
    mix[TUMOR] += share
    whole, single = [], []
    for present, rate in zip(_SIM["present"], _SIM["rates"]):
        n, m = present.shape
        cell = rng.choice(len(mix), size=n, p=mix)
        r = rate[cell]
        on = rng.random((n, m)) < r
        mean_rate = (r * present).sum(1) / np.maximum(present.sum(1), 1)
        locked = rng.random(n) < rho
        on[locked] = (rng.random(n) < mean_rate)[locked, None]
        x = np.where(present, np.where(on, 2, 1), 0)
        whole.append(x)
        rows, cols = np.nonzero(x)
        s = np.zeros((len(rows), m), dtype=x.dtype)
        s[np.arange(len(rows)), cols] = x[rows, cols]
        single.append(s)
    return whole, single


def sim_job(job):
    rho, share, rep, stack = job
    rng = np.random.default_rng([int(rho * 100), int(share * 1000), rep, stack])
    whole, single = [], []
    for _ in range(stack):
        w_, s_ = make_file(rho, share, rng)
        whole += w_
        single += s_
    w = LinkedISH(whole, _SIM["beta"] * stack)
    a_lnk, r_hat, _, conv, _ = w.fit(ish_start(w.t), 0.5, True, LINKED_ITERATIONS, ISH_STOP)
    if stack > 1:
        return (rho, share, stack), (r_hat, conv, a_lnk[EPI_MASK].sum(), np.nan, np.nan)
    s = LinkedISH(single, _SIM["beta"])
    assert w.n_beads.sum() == s.n_beads.sum() and not s.multi.any()
    a_ish = w.fit(ish_start(w.t), 0.0, False, ISH_ITERATIONS, ISH_STOP)[0]
    a_one = s.fit(ish_start(s.t), 0.0, False, ISH_ITERATIONS, ISH_STOP)[0]
    return (rho, share, stack), (r_hat, conv, a_lnk[EPI_MASK].sum(), a_ish[EPI_MASK].sum(), a_one[EPI_MASK].sum())


def made_up():
    sim_setup()
    jobs = [(r, s, k, 10) for r in RHOS for s in SHARES_STACKED for k in range(STACKS)]
    jobs += [(r, s, k, 1) for r in RHOS for s in SHARES for k in range(REPS)]
    with multiprocessing.Pool(12, initializer=sim_setup) as pool:
        res = pool.map(sim_job, jobs)
    ok = all(v[1] for _, v in res)
    print(f"\n  6  made-up data from the linked model, string places of {_SIM['ctrl']}.")
    print("     epithelial share = sum of the 16 epithelial types; the truth is the tumor share (all Colon-Ep).")
    print(f"\n  6a ten files stacked into one, {STACKS} stacks per row")
    print("true rho\ttrue tumor share\tfitted rho mean\tfitted rho sd\tlinked whole: epithelial mean")
    for r in RHOS:
        for s in SHARES_STACKED:
            v = np.array([x[1] for x in res if x[0] == (r, s, 10)], dtype=float)
            print(f"{r}\t{s}\t{v[:, 0].mean():.3f}\t{v[:, 0].std(ddof=1):.3f}\t{v[:, 2].mean():.4f}")
            ok &= abs(v[:, 0].mean() - r) <= 0.01
    print(f"\n  6b one real-sized file, {REPS} files per row")
    print("true rho\ttrue tumor share\tfitted rho mean\tfitted rho sd\tlinked whole: epithelial mean\tabs error\t"
          "CelFiE-ISH whole: epithelial mean\tabs error\tCelFiE-ISH single: epithelial mean\tabs error")
    for r in RHOS:
        for s in SHARES:
            v = np.array([x[1] for x in res if x[0] == (r, s, 1)], dtype=float)
            err = np.abs(v[:, 2:] - s).mean(0)
            print(f"{r}\t{s}\t{v[:, 0].mean():.3f}\t{v[:, 0].std(ddof=1):.3f}\t{v[:, 2].mean():.4f}\t{err[0]:.4f}\t"
                  f"{v[:, 3].mean():.4f}\t{err[1]:.4f}\t{v[:, 4].mean():.4f}\t{err[2]:.4f}")
            ok &= abs(v[:, 0].mean() - r) <= 0.02
    print(f"     6a fitted rho within 0.01 of the true rho in every row, 6b within 0.02, all fits converged: "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def main():
    os.makedirs("results/linked", exist_ok=True)
    print("EXPLORATORY: checks of the linked model's code (step_scripts/linked/check_model.py)")
    passed = []
    if sys.argv[1:] != ["made_up"]:
        with multiprocessing.Pool(6) as pool:
            for text, ok in pool.map(real_file_checks, six_files()):
                print(text)
                passed.append(ok)
    passed.append(made_up())
    print("\nALL CHECKS PASS" if all(passed) else "\nSOME CHECK FAILED")


if __name__ == "__main__":
    main()
