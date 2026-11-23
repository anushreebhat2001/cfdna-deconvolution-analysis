import sys

import numpy as np

import deconvolution_models.main as M
from deconvolution_models.celfie_ish import CelfieISH as Engine

out_unclipped, out_clipped = sys.argv[1], sys.argv[2]
np.random.seed(0)


def deconvolute(self):
    r = Engine(self.matrices, self.atlas_matrices, origins=self.origins,
               num_iterations=self.config["num_iterations"],
               convergence_criteria=self.config["stop_criterion"])
    self.alpha, self.i = r.two_step()
    self.alpha_clipped = r.get_clipped_alpha(0.05)


def write_output(self):
    np.savetxt(out_unclipped, (self.alpha * 100).reshape(1, -1), fmt="%.4f", delimiter="\t")
    np.savetxt(out_clipped, (self.alpha_clipped * 100).reshape(1, -1), fmt="%.4f", delimiter="\t")


M.CelfieISH.deconvolute = deconvolute
M.CelfieISH.write_output = write_output
M.main(["--model", "celfie-ish"] + sys.argv[3:] + ["--outfile", "/dev/null"], standalone_mode=False)
