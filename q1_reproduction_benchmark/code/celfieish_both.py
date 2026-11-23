"""usage: celfieish_both.py OUT_NOCLIP OUT_CLIP SEED <normal deconvolution arguments...>"""
import sys, numpy as np
out_noclip, out_clip, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
np.random.seed(seed)
import deconvolution_models.main as M
from deconvolution_models.celfie_ish import CelfieISH as Engine

def deconvolute(self):
    r = Engine(self.matrices, self.atlas_matrices, origins=self.origins,
               num_iterations=self.config["num_iterations"], convergence_criteria=self.config["stop_criterion"])
    self.alpha, self.i = r.two_step()
    self.alpha_clip = r.get_clipped_alpha(0.05)

def write_output(self):
    np.savetxt(out_noclip, (self.alpha * 100).reshape(1, -1), fmt="%.4f", delimiter="\t")
    np.savetxt(out_clip, (self.alpha_clip * 100).reshape(1, -1), fmt="%.4f", delimiter="\t")

M.CelfieISH.deconvolute = deconvolute
M.CelfieISH.write_output = write_output
M.main(sys.argv[4:], standalone_mode=False)
