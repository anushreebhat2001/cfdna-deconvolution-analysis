import argparse
import sys

import numpy as np

import deconvolution_models.main as M
from linked_model import LinkedISH

ISH_ITERATIONS, ISH_STOP = 10000, 0.0000001
LINKED_ITERATIONS = 200000


def load(pat_gz, atlas, regions):
    got = {}

    def deconvolute(self):
        got.update(matrices=self.matrices, atlas=self.atlas_matrices, origins=self.origins)

    old = M.CelfieISH.deconvolute, M.CelfieISH.write_output
    M.CelfieISH.deconvolute, M.CelfieISH.write_output = deconvolute, lambda self: None
    try:
        M.main(["--model", "celfie-ish", "-a", atlas, "-m", pat_gz, "--cpg_coordinates",
                "data/ref/hg38_pat_cpg.bed.gz", "-i", regions, "-b", "--epiformat", "pat",
                "--num_iterations", str(ISH_ITERATIONS), "--stop_criterion", str(ISH_STOP),
                "--random_restarts", "1", "--outfile", "/dev/null"], standalone_mode=False)
    finally:
        M.CelfieISH.deconvolute, M.CelfieISH.write_output = old
    return got["matrices"], got["atlas"], got["origins"]


def ish_start(n_types):
    np.random.seed(0)
    a = np.random.uniform(size=n_types)
    return a / a.sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stored")
    ap.add_argument("out")
    ap.add_argument("atlas")
    ap.add_argument("regions")
    ap.add_argument("--rho", type=float)
    args = ap.parse_args()

    eng = LinkedISH(*load(args.stored + ".pat.gz", args.atlas, args.regions))
    start = ish_start(eng.t)

    a0, _, it0, _, _ = eng.fit(start, 0.0, False, ISH_ITERATIONS, ISH_STOP)
    diff = {}
    for clip, ours in (("unclipped", a0), ("clipped", eng.clipped(a0, 0.0))):
        stored = np.loadtxt(f"{args.stored}.ish_{clip}.txt").ravel()
        diff[clip] = float(np.abs(np.round(ours * 100, 4) - stored).max())

    if args.rho is None:
        alpha, rho, it, ok, _ = eng.fit(start, 0.5, True, LINKED_ITERATIONS, ISH_STOP)
    else:
        alpha, rho, it, ok, _ = eng.fit(start, args.rho, False, LINKED_ITERATIONS, ISH_STOP)

    c0, c1 = eng.all_one_state_counts(a0, 0.0), eng.all_one_state_counts(alpha, rho)
    info = [("rho", f"{rho:.6f}"), ("rho_fitted", "no" if args.rho is not None else "yes"),
            ("iterations", it), ("converged", "yes" if ok else "no"),
            ("strings", len(eng.n_beads)), ("strings_2plus_beads", c0["multi"]),
            ("loglik_celfie_ish", f"{eng.loglik(a0, 0.0):.4f}"), ("loglik_linked", f"{eng.loglik(alpha, rho):.4f}"),
            ("celfie_ish_refit_iterations", it0),
            ("refit_max_diff_to_stored_unclipped", f"{diff['unclipped']:.4f}"),
            ("refit_max_diff_to_stored_clipped", f"{diff['clipped']:.4f}"),
            ("all_off_2plus_seen", c0["all_off"]), ("all_off_2plus_expected_celfie_ish", f"{c0['expected_all_off']:.2f}"),
            ("all_off_2plus_expected_linked", f"{c1['expected_all_off']:.2f}"),
            ("all_on_2plus_seen", c0["all_on"]), ("all_on_2plus_expected_celfie_ish", f"{c0['expected_all_on']:.2f}"),
            ("all_on_2plus_expected_linked", f"{c1['expected_all_on']:.2f}")]
    np.savetxt(args.out + ".ish_unclipped.txt", (alpha * 100).reshape(1, -1), fmt="%.4f", delimiter="\t")
    np.savetxt(args.out + ".ish_clipped.txt", (eng.clipped(alpha, rho) * 100).reshape(1, -1), fmt="%.4f", delimiter="\t")
    with open(args.out + ".info.tsv", "w") as fh:
        fh.writelines(f"{k}\t{v}\n" for k, v in info)


if __name__ == "__main__":
    sys.exit(main())
