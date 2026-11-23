"""Part A figures + method comparison table, from results/summary (run after collect_all.py).
Reads summary.txt, lod_bins.tsv, all_estimates.tsv; writes PNG + PDF figures and the table to SUMMARY_DIR/figures/.
usage: plot_part_a.py [SUMMARY_DIR]"""
import os, re, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

S = sys.argv[1] if len(sys.argv) > 1 else "/scratch/mxf4642/results/summary"
OUT = os.path.join(S, "figures")
os.makedirs(OUT, exist_ok=True)

# colours: chart chrome + categorical slots in fixed order (validated palette)
SURF, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
METHODS = [  # key, label, colour, marker  (colour follows the method in every figure)
    ("celfieish-clip",   "CelFiE-ISH, clipping",    "#2a78d6", "o"),
    ("celfieish-noclip", "CelFiE-ISH, no clipping", "#eb6834", "s"),
    ("celfeer",          "CelFEER",                 "#1baf7a", "^"),
    ("uxm",              "UXM",                     "#eda100", "D"),
    ("celfie",           "CelFiE",                  "#e87ba4", "v"),
]
LABEL = {k: l for k, l, _, _ in METHODS}
LABEL["celfie-orig"] = "CelFiE (original celfie.py)"
MARKERS = [("max", "Max markers (pan-epithelial, 8 groups)"), ("u250", "U250 markers (31 cell types)")]
SYNTH = ("uxm", "celfeer")          # on max markers these use synthetic megakaryocyte reads
THRESH = 0.95
FOOT_SYNTH = "* UXM and CelFEER on max markers use synthetic megakaryocyte reads (independent-sites approximation)."

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.titlecolor": INK, "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "axes.facecolor": SURF, "figure.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.7, "grid.linestyle": "-",
    "xtick.color": AXIS, "ytick.color": AXIS, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "xtick.major.size": 3, "ytick.major.size": 0, "legend.frameon": False, "legend.labelcolor": INK,
    "lines.linewidth": 1.6, "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round",
})


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.join(OUT, name + ".png/.pdf"))


def line(ax, x, y, col, mk, label=None, z=3):
    ax.plot(x, y, color=col, marker=mk, ms=5.5, mec=SURF, mew=1.2, label=label, zorder=z)


def thresh_handle():
    return Line2D([], [], color=MUTED, lw=0.9, label=f"Floor threshold (AUROC {THRESH})")


# ---------------------------------------------------------------- data
lod = pd.read_csv(f"{S}/lod_bins.tsv", sep="\t")
ref = lod[(lod.method == "celfieish-clip") & (lod.markers == "max")].sort_values("tf_lo")
BINLAB = [b.replace("%", "").replace("-", "\u2013") for b in ref["bin"]]
NPAT = ref["n_patients"].tolist()
X = np.arange(len(BINLAB))


def bins(key, ms):
    return lod[(lod.method == key) & (lod.markers == ms)].sort_values("tf_lo")


def bin_axis(ax, rot=40, fs=8, show=True):
    ax.set_xticks(X)
    ax.set_xticklabels(BINLAB if show else [], rotation=rot, ha="right" if rot else "center", fontsize=fs)
    ax.set_xlim(-0.4, len(X) - 0.6)


# ---------------------------------------------------------------- figure 1: detection-floor curves
fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
for ax, (ms, title) in zip(axes, MARKERS):
    for key, lab, col, mk in METHODS:
        d = bins(key, ms)
        if len(d):
            line(ax, X, d.auroc, col, mk, lab + (" *" if ms == "max" and key in SYNTH else ""))
    ax.axhline(THRESH, color=MUTED, lw=0.9, zorder=2)
    ax.set_title(title)
    bin_axis(ax)
    ax.set_xlabel("Tumor-fraction bin, % (12 virtual CRC samples per bin)")
axes[0].set_ylabel("AUROC vs 22 controls")
axes[0].set_ylim(0.35, 1.03)
h, l = axes[0].get_legend_handles_labels()
fig.legend(h + [thresh_handle()], l + [thresh_handle().get_label()], loc="lower center", ncol=3,
           bbox_to_anchor=(0.5, -0.17), handlelength=2.2, columnspacing=1.8)
fig.text(0.01, -0.21, FOOT_SYNTH, color=INK2, fontsize=8)
fig.suptitle("Detection floor by method, setup (a): AUROC per tumor-fraction bin", x=0.01, ha="left",
             fontweight="bold", color=INK, fontsize=11, y=1.02)
save(fig, "fig1_detection_floor")

# ---------------------------------------------------------------- figure 2: same, with bootstrap CIs (small multiples)
fig, axes = plt.subplots(2, len(METHODS), figsize=(15, 6.2), sharey=True, sharex=True)
for r, (ms, title) in enumerate(MARKERS):
    for c, (key, lab, col, mk) in enumerate(METHODS):
        ax = axes[r, c]
        d = bins(key, ms)
        if len(d):
            ax.fill_between(X, d.ci_lo, d.ci_hi, color=col, alpha=0.12, lw=0, zorder=1)
            line(ax, X, d.auroc, col, mk)
        ax.axhline(THRESH, color=MUTED, lw=0.9, zorder=2)
        if r == 0:
            ax.set_title(lab + (" *" if key in SYNTH else ""), fontsize=9.5)
        bin_axis(ax, rot=90, fs=7, show=(r == 1))
    axes[r, 0].set_ylabel(("Max markers" if ms == "max" else "U250 markers") + "\nAUROC (95% CI)")
axes[0, 0].set_ylim(0.0, 1.03)
fig.text(0.5, -0.03, "Tumor-fraction bin, %   (band = 95% bootstrap CI over patients and controls; "
         "bins hold 3\u20135 patients)", ha="center", color=INK2)
fig.text(0.01, -0.07, FOOT_SYNTH + f"   Gray line = floor threshold (AUROC {THRESH}).", color=INK2, fontsize=8)
fig.suptitle("Detection floor with 95% CIs, setup (a)", x=0.01, ha="left", fontweight="bold", color=INK,
             fontsize=11, y=1.01)
save(fig, "fig2_detection_floor_ci")

# ---------------------------------------------------------------- figure 3: dilution tracking
est = pd.read_csv(f"{S}/all_estimates.tsv", sep="\t")
crc = est[est.group == "crc"]
DILS = sorted([d for d in crc.dilution.unique() if d < 1], reverse=True)
XD = np.arange(len(DILS))
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
ymax = 1.5
for ax, (ms, title) in zip(axes, MARKERS):
    for key, lab, col, mk in METHODS:
        g = crc[(crc.method == key) & (crc.markers == ms)]
        if g.empty:
            continue
        base = g[g.dilution == 1.0].set_index("sample").epi
        base = base[base > 0]
        m = g[g.dilution < 1].join(base.rename("base"), on="sample", how="inner")
        m = m.assign(ratio=m.epi / (m.dilution * m.base))
        med = m.groupby("dilution").ratio.median().reindex(DILS)
        ymax = max(ymax, float(np.nanmax(med.values)))
        line(ax, XD, med.values, col, mk, lab + (" *" if ms == "max" and key in SYNTH else ""))
    ax.axhline(1.0, color=MUTED, lw=0.9, zorder=2)
    ax.set_title(title)
    ax.set_xticks(XD)
    ax.set_xticklabels([f"{100 * d:g}" for d in DILS])
    ax.set_xlabel("Share of the CRC sample in the mixture, %")
axes[0].set_ylim(0, min(ymax * 1.08, 4))
axes[0].set_ylabel("Measured / expected epithelial fraction\n(median over 20 CRC samples)")
h, l = axes[0].get_legend_handles_labels()
perfect = Line2D([], [], color=MUTED, lw=0.9, label="Perfect tracking (= 1)")
fig.legend(h + [perfect], l + [perfect.get_label()], loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.12),
           handlelength=2.2, columnspacing=1.8)
fig.text(0.01, -0.16, "Expected = dilution \u00d7 the undiluted sample's own estimate. " + FOOT_SYNTH,
         color=INK2, fontsize=8)
fig.suptitle("How well each method follows the two-fold dilution series, setup (a)", x=0.01, ha="left",
             fontweight="bold", color=INK, fontsize=11, y=1.02)
save(fig, "fig3_dilution_tracking")

# ---------------------------------------------------------------- figure 4: synthetic-read calibration (U250)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.0), sharey=True)
for ax, key in zip(axes, ("uxm", "celfeer")):
    for ms, lab, col, mk in (("u250", "Real Blood-Granul reference", "#2a78d6", "o"),
                             ("u250cal", "Synthetic Blood-Granul reference (independent sites)", "#eb6834", "s")):
        d = bins(key, ms)
        if len(d):
            line(ax, X, d.auroc, col, mk, lab)
    ax.axhline(THRESH, color=MUTED, lw=0.9, zorder=2)
    ax.set_title(f"{LABEL[key]}, U250 markers")
    bin_axis(ax)
    ax.set_xlabel("Tumor-fraction bin, %")
axes[0].set_ylabel("AUROC vs 22 controls")
axes[0].set_ylim(0.35, 1.03)
h, l = axes[0].get_legend_handles_labels()
fig.legend(h + [thresh_handle()], l + [thresh_handle().get_label()], loc="lower center", ncol=3,
           bbox_to_anchor=(0.5, -0.14), handlelength=2.2)
fig.suptitle("Calibration: effect of replacing real reads with independent-site synthetic reads",
             x=0.01, ha="left", fontweight="bold", color=INK, fontsize=11, y=1.03)
save(fig, "fig4_synthetic_calibration")

# ---------------------------------------------------------------- comparison table (parsed from summary.txt)
txt = open(f"{S}/summary.txt").read().split("\n")
R1 = re.compile(r"^\s+(\S+)\s+(\S+)\s+AUROC=([\d.]+) \[([\d.]+)-([\d.]+)\]\s+RMSE=([\d.]+)\s+R2=([\d.]+)\s+control median=([\d.]+)")
R2b = re.compile(r"^\s+(\S+)\s+(\S+)\s+sustained=(\S+)\s+first crossing=(\S+)")
R3 = re.compile(r"^\s+(\S+)\s+(\S+)\s+A=([\d.]+)\s+C=([\d.]+)\s+\(n=(\d+)\)")
R4 = re.compile(r"^\s+(\S+)\s+(\S+)\s+spearman\(depth, epi\)=([+-]?[\d.]+) \(p=([\d.]+)\)")
rows = {}
for t in txt:
    for rx, cols in ((R1, ["auroc", "ci_lo", "ci_hi", "rmse", "r2", "control_median"]),
                     (R2b, ["floor_sustained", "floor_first"]),
                     (R3, ["diag_A", "diag_C", "diag_n"]),
                     (R4, ["depth_rho", "depth_p"])):
        m = rx.match(t)
        if m:
            rows.setdefault((m.group(1), m.group(2)), {}).update(dict(zip(cols, m.groups()[2:])))
            break
PAPER = {("celfieish-clip", "u250"): "0.099", ("celfieish-noclip", "u250"): "0.126", ("celfeer", "u250"): "0.125",
         ("uxm", "u250"): "0.195", ("celfie", "u250"): "0.104", ("celfieish-clip", "max"): "~0.062",
         ("celfieish-noclip", "max"): "~0.066", ("celfie", "max"): "~0.072"}
NOTE = {("uxm", "max"): "synthetic megakaryocyte reads", ("celfeer", "max"): "synthetic megakaryocyte reads",
        ("uxm", "u250cal"): "calibration: synthetic Blood-Granul", ("celfeer", "u250cal"): "calibration: synthetic Blood-Granul",
        ("celfie-orig", "max"): "check: original celfie.py", ("celfie-orig", "u250"): "check: original celfie.py"}
order_m = [k for k, _, _, _ in METHODS] + ["celfie-orig"]
order_s = ["max", "u250", "u250cal"]
tab = pd.DataFrame([dict(method=LABEL.get(k[0], k[0]), markers=k[1], **v, paper_rmse=PAPER.get(k, ""),
                         note=NOTE.get(k, ""), _o=(order_s.index(k[1]) if k[1] in order_s else 9,
                                                   order_m.index(k[0]) if k[0] in order_m else 9))
                    for k, v in rows.items()]).sort_values("_o").drop(columns="_o")
tab.to_csv(os.path.join(OUT, "table_method_comparison.tsv"), sep="\t", index=False)
md = ["| Method | Markers | AUROC [95% CI] | RMSE (paper) | R\u00b2 | Control median | Sustained floor | "
      "First crossing | Diagonal RMSE (C) | Depth \u03c1 (p) | Note |", "|" + "---|" * 11]
for _, r in tab.iterrows():
    md.append(f"| {r.method} | {r.markers} | {r.auroc} [{r.ci_lo}\u2013{r.ci_hi}] | {r.rmse}"
              f"{' (' + r.paper_rmse + ')' if r.paper_rmse else ''} | {r.r2} | {r.control_median} | "
              f"{r.floor_sustained} | {r.floor_first} | {r.diag_C} | {r.depth_rho} ({r.depth_p}) | {r.note} |")
open(os.path.join(OUT, "table_method_comparison.md"), "w").write("\n".join(md) + "\n")
print("wrote", os.path.join(OUT, "table_method_comparison.tsv/.md"), f"({len(tab)} rows)")
print("\n".join(md))
