import collections
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
sys.path.insert(0, os.path.join(HERE, "..", "factors"))
import score_pair
import stage1_checks as s1

L = "results/linked"
META, auroc = s1.META, s1.auroc
CTRL = [m for m in META if not m["cancer"]]
CANC = [m for m in META if m["cancer"]]
EPI_MASK = s1.MASK["epithelial16"]
SETS = [("A1 all strings", "step2", "whole"), ("A2 strings with 2+ beads", "step2_2plus", "whole"),
        ("A3 4+ bead strings, whole", "step4", "whole"), ("A3 4+ bead strings, 3-bead pieces", "step4", "size3"),
        ("A3 4+ bead strings, 2-bead pieces", "step4", "size2"), ("A4 2-bead strings", "step3", "b2_whole"),
        ("A4 3-bead strings", "step3", "b3_whole"), ("A5 bead-matched buckets", "step3_matched", "d*_b*_whole"),
        ("B3 HU012 background", "b3/HU012", "whole"), ("B3 HU083 background", "b3/HU083", "whole"),
        ("D1 simulated", "d1", "whole"), ("breast", "breast", "whole"), ("pan-epithelial (step 7)", "step7", "whole")]
SETS += [(f"T10 {c}", f"factors/t10_{c}", "whole") for c in
         ("scattered_1", "clumped_1", "scattered_3", "clumped_3", "symmetric_5", "symmetric_10")]


def info(path):
    return dict(line.rstrip("\n").split("\t") for line in open(path))


def epi(path, mask=EPI_MASK):
    v = np.loadtxt(path).ravel()
    return float(v[mask].sum() / v.sum())


def arms(folder, i, clip="unclipped", single="cut", whole="whole"):
    return np.array([epi(f"{L}/{folder}/{i}/{whole}_linked.ish_{clip}.txt"),
                     epi(f"results/{folder}/{i}/{whole}.ish_{clip}.txt"),
                     epi(f"results/{folder}/{i}/{single}.ish_{clip}.txt")])


def health():
    print("1. The run: every linked fit also refitted its file with rho = 0 and compared with the stored CelFiE-ISH answer")
    print("set\tfits\tfailed\tnot converged\tlargest difference of the rho = 0 refit to stored (percent points)\t"
          "fits where the linked likelihood is below CelFiE-ISH's\tmedian EM rounds\tmost EM rounds")
    for name, folder, variant in SETS:
        files = [f for f in glob.glob(f"{L}/{folder}/*/{variant}_linked.info.tsv") if not os.path.islink(os.path.dirname(f))]
        inf = [info(f) for f in files]
        failed = len(glob.glob(f"{L}/{folder}/*/{variant}_linked.failed"))
        if not inf:
            print(f"{name}\t0\t{failed}")
            continue
        diff = max(max(float(x["refit_max_diff_to_stored_unclipped"]), float(x["refit_max_diff_to_stored_clipped"])) for x in inf)
        worse = sum(float(x["loglik_linked"]) < float(x["loglik_celfie_ish"]) - 1e-4 for x in inf)
        its = [int(x["iterations"]) for x in inf]
        print(f"{name}\t{len(inf)}\t{failed}\t{sum(x['converged'] != 'yes' for x in inf)}\t{diff:.4f}\t{worse}\t"
              f"{int(np.median(its))}\t{max(its)}")
    same = glob.glob(f"{L}/*/*/whole_samerho.info.tsv") + glob.glob(f"{L}/b3/*/*__d*/whole_samerho.info.tsv")
    print(f"same-rho fits\t{len(same)}\t{len(glob.glob(f'{L}/**/*_samerho.failed', recursive=True))}\t"
          f"{sum(info(f)['converged'] != 'yes' for f in same)}")


def rho_table():
    print("\n\n2. Fitted rho (share of 2+ bead strings in the locked part)")
    print("files\tn\tmedian\tlowest\thighest")

    def row(name, paths):
        r = [float(info(p)["rho"]) for p in paths]
        print(f"{name}\t{len(r)}\t{np.median(r):.3f}\t{min(r):.3f}\t{max(r):.3f}")

    f = L + "/step2/{}/whole_linked.info.tsv"
    for b in ("HU012", "HU083"):
        row(f"A1 healthy controls {b}", [f.format(m["id"]) for m in CTRL if m["batch"] == b])
    for b in sorted({m["batch"] for m in CANC}):
        row(f"A1 cancer samples run {b}", [f.format(m["id"]) for m in CANC if m["batch"] == b])
    row("A1 dilutions", glob.glob(f"{L}/step2/*__d*/whole_linked.info.tsv"))
    row("breast samples", glob.glob(f"{L}/breast/*/whole_linked.info.tsv"))
    row("A3 4+ bead strings, 42 study samples", [f"{L}/step4/{m['id']}/whole_linked.info.tsv" for m in META])
    row("A4 2-bead strings, 42 study samples", [f"{L}/step3/{m['id']}/b2_whole_linked.info.tsv" for m in META])
    row("A4 3-bead strings, 42 study samples", [f"{L}/step3/{m['id']}/b3_whole_linked.info.tsv" for m in META])
    row("pan-epithelial, 42 study samples", [f"{L}/step7/{m['id']}/whole_linked.info.tsv" for m in META])
    row("D1 simulated (beads unrelated), all 142", glob.glob(f"{L}/d1/*/whole_linked.info.tsv"))
    for c in ("scattered_1", "clumped_1", "scattered_3", "clumped_3", "symmetric_5", "symmetric_10"):
        row(f"T10 {c}", glob.glob(f"{L}/factors/t10_{c}/*/whole_linked.info.tsv"))
    print(f"rho used in the same-rho check: {open(L + '/same_rho.txt').read().strip()}")


def all_one_state():
    print("\n\n3. Strings with 2+ beads that are all-off or all-on: seen, and expected by each fitted model (A1, summed over files)")
    print("files\tstrings with 2+ beads\tall-off seen\texpected CelFiE-ISH\texpected linked\tall-on seen\texpected CelFiE-ISH\texpected linked")
    groups = [(f"healthy controls {b}", [m["id"] for m in CTRL if m["batch"] == b]) for b in ("HU012", "HU083")]
    groups += [("cancer samples", [m["id"] for m in CANC]),
               ("D1 simulated controls", None)]
    keys = ("strings_2plus_beads", "all_off_2plus_seen", "all_off_2plus_expected_celfie_ish", "all_off_2plus_expected_linked",
            "all_on_2plus_seen", "all_on_2plus_expected_celfie_ish", "all_on_2plus_expected_linked")
    for name, ids in groups:
        folder = "step2" if ids else "d1"
        inf = [info(f"{L}/{folder}/{i}/whole_linked.info.tsv") for i in (ids or [m["id"] for m in CTRL])]
        t = [sum(float(x[k]) for x in inf) for k in keys]
        print(f"{name}\t" + "\t".join(f"{v:.0f}" for v in t))


def healthy(folder="step2", label="A1 all strings", single="cut", whole="whole"):
    rng = np.random.default_rng(0)
    print(f"\n{label}")
    print("clip\tmedian healthy epithelial share HU012: linked\tISH whole\tsingle\tHU083: linked\tISH whole\tsingle\t"
          "healthy HU012 against healthy HU083 (AUROC): linked\tISH whole\tsingle\tlinked minus ISH whole\tci_lo\tci_hi")
    for clip in ("unclipped", "clipped"):
        sc = {b: np.array([arms(folder, m["id"], clip, single, whole) for m in CTRL if m["batch"] == b]) for b in ("HU012", "HU083")}
        med = [np.median(sc[b][:, k]) for b in ("HU012", "HU083") for k in range(3)]
        au = [auroc(sc["HU012"][:, k], sc["HU083"][:, k]) for k in range(3)]
        bs = []
        for _ in range(s1.N_BOOT):
            p = sc["HU012"][rng.integers(0, len(sc["HU012"]), len(sc["HU012"]))]
            n = sc["HU083"][rng.integers(0, len(sc["HU083"]), len(sc["HU083"]))]
            bs.append(auroc(p[:, 0], n[:, 0]) - auroc(p[:, 1], n[:, 1]))
        lo, hi = np.percentile(bs, [2.5, 97.5])
        print(f"{clip}\t" + "\t".join(f"{v:.3f}" for v in med + au) + f"\t{au[0] - au[1]:+.3f}\t{lo:+.3f}\t{hi:+.3f}")


def summary(label):
    out, clip, cols = {}, None, None
    for line in open(f"{L}/scores/{label}_summary.tsv"):
        f = line.rstrip("\n").split("\t")
        if line.startswith("==="):
            clip = "clipped" if "ISH clipped" in line else "unclipped"
        elif f[0] == "subset":
            cols = f
        elif cols and len(f) == len(cols):
            out[(clip, f[0], f[1])] = dict(zip(cols, f))
    return out


def detection(title, vs_single, vs_ish, subsets, controls=("all", "HU012", "HU083")):
    a, b = summary(vs_single), summary(vs_ish)
    print(f"\n{title}")
    print("clip\tcancer files\tcontrols\tn files\tAUROC linked\tISH whole\tsingle\tlinked minus single\tci_lo\tci_hi\tmoved\t"
          "linked minus ISH whole\tmoved\terror linked\tISH whole\tsingle\tsingle minus linked\tmoved")
    for clip in ("unclipped", "clipped"):
        for sub in subsets:
            for c in controls:
                if (clip, sub, c) not in a:
                    continue
                r, q = a[(clip, sub, c)], b[(clip, sub, c)]
                print(f"{clip}\t{sub}\t{c}\t{r['cancer_files']}\t{r['auroc_A']}\t{q['auroc_B']}\t{r['auroc_B']}\t{r['auroc_gain']}\t"
                      f"{r['auc_ci_lo']}\t{r['auc_ci_hi']}\t{r['auroc_moved']}\t{q['auroc_gain']}\t{q['auroc_moved']}\t"
                      f"{r['error_A']}\t{q['error_B']}\t{r['error_B']}\t{r['error_gain']}\t{r['error_moved']}")


def low6_and_breast():
    print("\n\n6. Low-tumor samples and the breast samples (A1, undiluted)")
    print("clip\tpositives\tcontrols\tAUROC linked\tsingle\tgain\tci_lo\tci_hi\tmoved\tAUROC ISH whole")
    for clip in ("unclipped", "clipped"):
        sc = {m["id"]: arms("step2", m["id"], clip) for m in META}
        br = {i: arms("breast", i, clip) for i in s1.BREAST}
        for gname, groups in (("low6 (tumor share under 20%)", [[sc[m["id"]] for m in CANC if m["tf"] < s1.LOW and m["patient"] == p]
                                                                 for p in sorted({m["patient"] for m in CANC})]),
                              ("breast4", [[v] for v in br.values()])):
            groups = [np.array(g) for g in groups if g]
            for k in ("all", "HU012", "HU083"):
                neg = np.array([sc[m["id"]] for m in CTRL if k in ("all", m["batch"])])
                rng = np.random.default_rng(0)
                r = s1.boot_gain([g[:, [0, 2]] for g in groups], neg[:, [0, 2]], rng)
                pos = np.concatenate(groups)
                print(f"{clip}\t{gname}\t{k}\t{s1.fmt(r)}\t{auroc(pos[:, 1], neg[:, 1]):.3f}")
    print("\nbreast samples, unclipped epithelial share and the number of the 22 controls each one is above")
    print("sample\tlinked\tabove\tISH whole\tabove\tsingle\tabove")
    sc = np.array([arms("step2", m["id"]) for m in CTRL])
    for i in s1.BREAST:
        v = arms("breast", i)
        print(f"{i}\t" + "\t".join(f"{v[k]:.3f}\t{int((v[k] > sc[:, k]).sum())}" for k in range(3)))


def t10():
    print("\n\n7. T10: simulated healthy controls (true epithelial share 0) with one switch on. Epithelial score = false signal.")
    ids = [m["id"] for m in CTRL]
    for clip in ("unclipped", "clipped"):
        print(f"\n{clip}\ncondition\tmedian fitted rho\tmedian false signal: linked\tISH whole\tsingle\tmean linked minus single\tci_lo\tci_hi\tmoved")
        rng = np.random.default_rng(0)
        for name, folder in [("none (D1 as stored)", "d1")] + [(c, f"factors/t10_{c}") for c in
                             ("scattered_1", "clumped_1", "scattered_3", "clumped_3", "symmetric_5", "symmetric_10")]:
            a = np.array([arms(folder, i, clip) for i in ids])
            rho = np.median([float(info(f"{L}/{folder}/{i}/whole_linked.info.tsv")["rho"]) for i in ids])
            d = a[:, 0] - a[:, 2]
            bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(s1.N_BOOT)]
            lo, hi = np.percentile(bs, [2.5, 97.5])
            print(f"{name}\t{rho:.3f}\t" + "\t".join(f"{np.median(a[:, k]):.3f}" for k in range(3))
                  + f"\t{d.mean():+.3f}\t{lo:+.3f}\t{hi:+.3f}\t{'yes' if lo > 0 or hi < 0 else 'no'}")


def a5():
    print("\n\n8. CpGs per read with equal bead totals (A5, 10 draws): linked whole against single, undiluted samples and all files")
    print("clip\tbeads per string\tcancer files\tcontrols\tAUROC linked\tsingle\tgain\tci_lo\tci_hi\tmoved\terror linked\tsingle\tgain\tmoved\t"
          "for comparison CelFiE-ISH whole against single: AUROC gain\tmoved")
    for clip in ("unclipped", "clipped"):
        for k in ("2", "3", "4+"):
            a = summary(f"A5_matched_bucket{k}_linked_vs_single")
            old = collections.defaultdict(dict)
            c2 = None
            for line in open(f"results/scores/A5_matched_bucket{k}_summary.tsv"):
                f = line.rstrip("\n").split("\t")
                if line.startswith("==="):
                    c2 = "clipped" if "ISH clipped" in line else "unclipped"
                elif len(f) > 16 and f[0] != "subset":
                    old[(c2, f[0], f[1])] = (f[13], f[16])
            for sub in ("dilution_1", "all"):
                for c in ("all", "HU012", "HU083"):
                    r = a[(clip, sub, c)]
                    print(f"{clip}\t{k}\t{sub}\t{c}\t{r['auroc_A']}\t{r['auroc_B']}\t{r['auroc_gain']}\t{r['auc_ci_lo']}\t{r['auc_ci_hi']}\t"
                          f"{r['auroc_moved']}\t{r['error_A']}\t{r['error_B']}\t{r['error_gain']}\t{r['error_moved']}\t"
                          + "\t".join(old[(clip, sub, c)]))


def main():
    print("EXPLORATORY: the linked model (CelFiE-ISH plus one togetherness number, rho). Plan: step_scripts/linked/PLAN.md")
    print("Code checks: results/linked/model_checks.txt\n")
    health()
    rho_table()
    all_one_state()
    print("\n\n4. Healthy controls: false epithelial signal, and how well it separates the two healthy runs (L3)")
    healthy()
    healthy("step2_2plus", "A2 strings with 2+ beads")
    healthy("step4", "A3 4+ bead strings", single="size1")
    low = ("all", "dilution_1", "tumor_below_1pct", "tumor_1_to_5pct", "tumor_5_to_20pct", "tumor_20pct_up")
    print("\n\n5. Detection and tumor amount, three arms (L4, L5). Gain intervals: patients and controls redrawn, 2000 rounds.")
    detection("A1 all strings", "A1_linked_vs_single", "A1_linked_vs_ish", low)
    detection("A1 all strings, SAME rho for every file", "A1_samerho_vs_single", "A1_linked_vs_ish", low)
    detection("A2 strings with 2+ beads", "A2_linked_vs_single", "A2_linked_vs_ish", low)
    detection("A3 4+ bead strings", "A3_linked_whole_vs_single", "A3_linked_whole_vs_ish", low)
    for b in ("HU012", "HU083"):
        detection(f"B3 dilutions into {b} background (matched row: controls {b})", f"B3_background_{b}_linked_vs_single",
                  f"B3_background_{b}_linked_vs_ish", low)
        detection(f"B3 dilutions into {b} background, SAME rho for every file", f"B3_background_{b}_samerho_vs_single",
                  f"B3_background_{b}_linked_vs_ish", low)
    detection("D1 simulated, beads unrelated (L1)", "D1_linked_vs_single", "D1_linked_vs_ish", low)
    detection("Pan-epithelial markers (step 7)", "PanEpi_linked_vs_single", "PanEpi_linked_vs_ish", low)
    low6_and_breast()
    t10()
    a5()


if __name__ == "__main__":
    main()
