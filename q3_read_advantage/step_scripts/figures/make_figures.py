import os
import sys
from decimal import ROUND_HALF_UP, Decimal

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("common", "checks", "factors"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import four_checks
import score_pair
import stage1_checks as s1
import stage2_readout as s2

OUT = "writeup/figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
WHOLE, SINGLE = "#2a78d6", "#eb6834"
plt.rcParams.update({"figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
                     "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.edgecolor": GRID, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
                     "axes.titlelocation": "left", "font.family": "DejaVu Sans"})


def clean(ax, grid="y"):
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.grid(axis=grid, color=GRID, linewidth=1)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)


def und_row(label, controls="all"):
    r = next(r for r in four_checks.summary(label)["unclipped"][0]
             if r["subset"] == "dilution_1" and r["controls_used"] == controls)
    return float(r["auroc_gain"]), float(r["auc_ci_lo"]), float(r["auc_ci_hi"])


def fig1():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.9), sharey=True, gridspec_kw={"wspace": 0.08})
    panels = [("Reads grouped by CpGs per read",
               "CpGs per read (inside a marker)", ["2", "3", "4 or more"],
               [("every read of the group", ["A4_bucket2", "A4_bucket3", "A3_ladder_whole_vs_single"], "o", -0.15),
                ("groups thinned to equal CpG totals (average of 10 draws)",
                 ["A5_matched_bucket2", "A5_matched_bucket3", "A5_matched_bucket4+"], "s", 0.15)]),
              ("Reads with 4+ CpGs, cut into pieces",
               "CpGs kept together", ["2", "3", "whole read"],
               [("every read of the group", ["A3_ladder_2bead_vs_single", "A3_ladder_3bead_vs_single",
                                               "A3_ladder_whole_vs_single"], "o", 0.0)])]
    handles = {}
    for ax, (title, xlabel, ticks, series) in zip(axes, panels):
        clean(ax)
        ax.axhline(0, color=INK2, linewidth=1)
        for name, labels, marker, dx in series:
            x = np.arange(len(labels)) + dx
            g = np.array([und_row(lab) for lab in labels])
            ax.errorbar(x, g[:, 0], yerr=[g[:, 0] - g[:, 1], g[:, 2] - g[:, 0]], fmt="none", ecolor=WHOLE, elinewidth=2,
                        alpha=0.45 if marker == "s" else 1)
            handles[name], = ax.plot(x, g[:, 0], marker, color=WHOLE, markersize=9 if marker == "o" else 7,
                                     markeredgecolor=SURFACE if marker == "o" else WHOLE, markeredgewidth=2,
                                     markerfacecolor=WHOLE if marker == "o" else SURFACE, linestyle="none")
            for xi, (v, lo, hi) in zip(x, g):
                side = 1 if marker == "s" or dx == 0 else -1
                ax.annotate(f"{v:+.2f}", (xi, v), xytext=(9 * side, 0), textcoords="offset points",
                            ha="left" if side > 0 else "right", va="center", color=INK, fontsize=9)
            for b, mk in (("HU012", "_"), ("HU083", "|")):
                handles[f"point value against {b} controls only"], = ax.plot(
                    x, [und_row(lab, b)[0] for lab in labels], mk, color=INK2, markersize=8, markeredgewidth=1.5, linestyle="none")
        ax.set_xticks(range(len(ticks)), ticks)
        ax.set_xlim(-0.5, len(ticks) - 0.5)
        ax.set_xlabel(xlabel)
        ax.set_title(title)
    axes[0].set_ylabel("gain in detection score (AUROC)\nwhole reads minus single CpGs")
    axes[0].set_ylim(-0.05, 0.32)
    fig.legend(list(handles.values()), list(handles), frameon=False, fontsize=8.5, loc="lower left", ncols=2,
               bbox_to_anchor=(0.065, 0.085), labelcolor=INK2, columnspacing=2.5)
    fig.suptitle("Planned result: gain from reading whole reads, by CpGs per read", x=0.065, ha="left", fontsize=13, fontweight="bold")
    fig.text(0.065, 0.008, "20 undiluted cancer samples against 22 healthy controls, U250 markers, unclipped. Vertical bars: 95% bootstrap "
             "interval (patients and controls redrawn).\nAbove zero = whole reads better. A gain here cannot be credited to tumor: see the "
             "healthy-against-healthy comparison in the next figure.", fontsize=8.5, color=INK2, va="bottom")
    fig.subplots_adjust(left=0.085, right=0.985, top=0.87, bottom=0.31)
    fig.savefig(f"{OUT}/fig1_gain_by_beads.png", dpi=200)
    plt.close(fig)


def detection_panel(ax, title, groups):
    clean(ax)
    w = 0.27
    for k, (label, a, b, gain) in enumerate(groups):
        if a is None:
            ax.text(k, 0.25, "not fitted\nyet", ha="center", color=INK2, fontsize=9)
            continue
        ax.bar(k - w / 2 - 0.015, a, w, color=WHOLE, label="whole reads" if k == 0 else None)
        ax.bar(k + w / 2 + 0.015, b, w, color=SINGLE, label="single CpGs" if k == 0 else None)
        ax.text(k - w / 2 - 0.015, a + 0.012, f"{a:.2f}", ha="center", fontsize=9, color=INK)
        ax.text(k + w / 2 + 0.015, b + 0.012, f"{b:.2f}", ha="center", fontsize=9, color=INK)
        if gain:
            g = [float(Decimal(repr(float(v))).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP)) for v in gain]
            ax.text(k, 1.105, f"gain {g[0]:+.2f}\n({g[1]:+.2f} to {g[2]:+.2f})", ha="center", fontsize=8, color=INK2)
    ax.axhline(0.5, color=INK2, linewidth=1)
    ax.set_xlim(-0.5, len(groups) - 0.1)
    ax.text(len(groups) - 0.12, 0.508, "chance", ha="right", va="bottom", fontsize=8, color=INK2)
    ax.set_xticks(range(len(groups)), [g[0] for g in groups], fontsize=8.5)
    ax.set_ylim(0, 1.2)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_title(title, pad=14)


def auc_rows(sc, breast):
    rng = np.random.default_rng(0)
    neg = np.array([sc[m["id"]] for m in s2.CTRL])
    out = []
    for keep in (lambda m: True, lambda m: m["tf"] < s1.LOW):
        by = {}
        for m in s2.CANC:
            if keep(m):
                by.setdefault(m["patient"], []).append(sc[m["id"]])
        out.append(s1.boot_gain([np.array(v) for v in by.values()], neg, rng))
    out.append(s1.boot_gain([sc[i][None, :] for i in s1.BREAST], neg, rng) if breast else None)
    h = {b: np.array([sc[m["id"]] for m in s2.CTRL if m["batch"] == b]) for b in ("HU012", "HU083")}
    out.append(s1.boot_gain([v[None, :] for v in h["HU012"]], h["HU083"], rng))
    return out


def fig2():
    ids = [m["id"] for m in s1.META]
    names = ["cancer, all 20\n(planned)", "cancer, tumor\nunder 20% (6)", "breast, almost\nno tumor (4)", "healthy HU012\nagainst HU083"]
    u = s2.scores(lambda i: "results/breast" if i in s1.BREAST else "results/step2", s1.MASK["epithelial16"], ids + s1.BREAST)
    mnames = [h[:-5] for h in open(score_pair.ATLAS["max"]).readline().rstrip("\n").split("\t") if h.endswith("_METH")]
    epi = np.array([n == "Epithelial" for n in mnames])
    has_breast = all(os.path.exists(f"results/breast_max/{i}/{v}.ish_unclipped.txt") for i in s1.BREAST for v in ("whole", "cut"))
    p = s2.scores(lambda i: "results/breast_max" if i in s1.BREAST else "results/step7", epi, ids + (s1.BREAST if has_breast else []))
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.0), sharey=True, gridspec_kw={"wspace": 0.06})
    for ax, title, sc, br in ((axes[0], "U250 markers (31 cell types)", u, True),
                              (axes[1], "Pan-epithelial markers (8 cell types)  ·  exploratory", p, has_breast)):
        rows = auc_rows(sc, br)
        if sc is u:
            rows[0] = (*rows[0][:2], *und_row("A1_all_strings"))
        detection_panel(ax, title, [(n, *(r[:2] if r else (None, None)), r[2:] if r else None) for n, r in zip(names, rows)])
    axes[0].set_ylabel("detection score (AUROC)")
    fig.legend(*axes[0].get_legend_handles_labels(), frameon=False, loc="lower left", fontsize=9, ncols=2,
               bbox_to_anchor=(0.055, 0.1), labelcolor=INK2)
    fig.suptitle("On U250, whole reads separate tumor-free groups almost as well as cancer", x=0.06, ha="left", fontsize=13, fontweight="bold")
    fig.text(0.06, 0.012, "First three groups: scored against the 22 healthy controls. Last group: the two healthy batches against each other "
             "(no tumor on either side). Unclipped.\nGain = whole minus single, with its 95% bootstrap interval. Only the first U250 group "
             "was planned; the others were chosen after seeing the results.", fontsize=8.5, color=INK2, va="bottom")
    fig.subplots_adjust(left=0.06, right=0.99, top=0.84, bottom=0.27)
    fig.savefig(f"{OUT}/fig2_detection.png", dpi=200)
    plt.close(fig)


def fig3():
    ids = [m["id"] for m in s2.CTRL]
    u = s2.scores(lambda i: "results/step2", s1.MASK["epithelial16"], ids)
    mnames = [h[:-5] for h in open(score_pair.ATLAS["max"]).readline().rstrip("\n").split("\t") if h.endswith("_METH")]
    p = s2.scores(lambda i: "results/step7", np.array([n == "Epithelial" for n in mnames]), ids)
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    clean(ax)
    rng = np.random.default_rng(1)
    groups = [("U250\nhealthy HU012", u, "HU012"), ("U250\nhealthy HU083", u, "HU083"),
              ("pan-epithelial\nhealthy HU012", p, "HU012"), ("pan-epithelial\nhealthy HU083", p, "HU083")]
    for k, (label, sc, batch) in enumerate(groups):
        v = 100 * np.array([sc[m["id"]] for m in s2.CTRL if m["batch"] == batch])
        for j, (color, dx, name) in enumerate(((WHOLE, -0.17, "whole reads"), (SINGLE, 0.17, "single CpGs"))):
            x = k + dx + rng.uniform(-0.06, 0.06, len(v))
            ax.plot(x, v[:, j], "o", color=color, markersize=6.5, markeredgecolor=SURFACE, markeredgewidth=1.2, alpha=0.9,
                    label=f"{name} (one dot per healthy person)" if k == 0 else None)
            med = np.median(v[:, j])
            ax.plot([k + dx - 0.13, k + dx + 0.13], [med, med], color=INK, linewidth=2, label="median" if k == 0 and j == 1 else None)
            ax.text(k + dx + (0.16 if j else -0.16), med, f"{med:.1f}", va="center", ha="left" if j else "right", fontsize=9, color=INK)
    ax.axvline(1.5, color=GRID, linewidth=1)
    ax.axhline(0, color=INK2, linewidth=1)
    ax.set_xticks(range(4), [g[0] for g in groups])
    ax.set_xlim(-0.6, 3.6)
    ax.set_ylabel("estimated epithelial share (%)\nhealthy blood should read about 0")
    ax.legend(frameon=False, fontsize=9, loc="upper right", labelcolor=INK2)
    ax.set_title("False epithelial signal in healthy controls, by marker set  ·  exploratory", pad=10)
    fig.text(0.075, 0.012, "22 healthy controls, all reads, unclipped. Removing the megakaryocyte column from the pan-epithelial atlas "
             "leaves its\nfloor unchanged (3% and 0%), so that column is not the reason for the difference.", fontsize=8.5, color=INK2, va="bottom")
    fig.subplots_adjust(left=0.105, right=0.985, top=0.9, bottom=0.215)
    fig.savefig(f"{OUT}/fig3_healthy_floor.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for name in sys.argv[1:] or ["fig1", "fig2", "fig3"]:
        globals()[name]()
        print("wrote", name)
