import numpy as np

from deconvolution_models.celfie_ish import CelfieISH


class LinkedISH(CelfieISH):
    def __init__(self, mixtures, beta, origins=None):
        super().__init__(mixtures, beta, origins=origins)
        self.L_ind = np.hstack(self.log_term1)
        lock, den, rates, n_on, n_off, all_on_ind, all_off_ind = [], [], [], [], [], [], []
        for x_m, x_u, b in zip(self.x_c_m, self.x_c_u, self.beta):
            x_m, x_u = x_m.astype(int), x_u.astype(int)
            x_v = (x_m + x_u).T
            d = (~np.isnan(b)).astype(int) @ x_v
            on, off = x_m.sum(1), x_u.sum(1)
            with np.errstate(invalid="ignore", divide="ignore"):
                mean_rate = (np.nan_to_num(b) @ x_v) / d
                lock.append(np.log(np.where(off == 0, mean_rate, np.where(on == 0, 1 - mean_rate, 0.0))))
                all_on_ind.append(np.nan_to_num(np.log(b)) @ x_v)
                all_off_ind.append(np.nan_to_num(np.log(1 - b)) @ x_v)
            den.append(d)
            rates.append(mean_rate)
            n_on.append(on)
            n_off.append(off)
        self.n_on, self.n_off = np.concatenate(n_on), np.concatenate(n_off)
        self.n_beads = self.n_on + self.n_off
        self.multi = self.n_beads >= 2
        self.den, self.mean_rate = np.hstack(den), np.hstack(rates)
        self.L_lock = np.hstack(lock)
        self.L_lock[self.den == 0] = 0.0
        self.L_lock[:, ~self.multi] = self.L_ind[:, ~self.multi]
        self.L_all_on_ind, self.L_all_off_ind = np.hstack(all_on_ind), np.hstack(all_off_ind)
        self.scale = np.maximum(self.L_ind, self.L_lock).max(axis=0)
        self.P_ind, self.P_lock = np.exp(self.L_ind - self.scale), np.exp(self.L_lock - self.scale)
        self.P_lock_multi = self.P_lock[:, self.multi]

    def joint(self, alpha, rho):
        j = alpha[:, None] * ((1 - rho) * self.P_ind + rho * self.P_lock)
        return j, j.sum(axis=0)

    def loglik(self, alpha, rho):
        return float((np.log(self.joint(alpha, rho)[1]) + self.scale).sum())

    def posterior(self, alpha, rho):
        j, norm = self.joint(alpha, rho)
        return j / norm

    def fit(self, alpha0, rho0, fit_rho, num_iterations, stop, rho_stop=1e-8, track=False):
        alpha, rho, trace = np.array(alpha0, dtype=float), float(rho0), []
        n_multi = int(self.multi.sum())
        for i in range(num_iterations):
            j, norm = self.joint(alpha, rho)
            if track:
                trace.append(float((np.log(norm) + self.scale).sum()))
            new_alpha = (j / norm).sum(axis=1)
            new_alpha /= new_alpha.sum()
            assert not np.isnan(new_alpha).any(), "alpha has NaN"
            new_rho = rho
            if fit_rho and n_multi:
                new_rho = float(rho * ((alpha @ self.P_lock_multi) / norm[self.multi]).sum() / n_multi)
            if i and np.mean(abs(new_alpha - alpha)) / np.mean(abs(alpha)) < stop and abs(new_rho - rho) < rho_stop:
                return alpha, rho, i, True, trace
            alpha, rho = new_alpha, new_rho
        return alpha, rho, num_iterations, False, trace

    def clipped(self, alpha, rho, clipping=0.05):
        z = self.posterior(alpha, rho)
        z[z < clipping] = 0
        out = z.sum(axis=1)
        return out / out.sum()

    def all_one_state_counts(self, alpha, rho):
        m = self.multi
        known = self.den[:, m] > 0
        w = alpha[:, None] * known
        w = w / w.sum(axis=0)
        rate = np.where(known, self.mean_rate[:, m], 0.5)
        p_off = (1 - rho) * np.exp(self.L_all_off_ind[:, m]) + rho * (1 - rate)
        p_on = (1 - rho) * np.exp(self.L_all_on_ind[:, m]) + rho * rate
        return dict(multi=int(m.sum()), all_off=int((self.n_on[m] == 0).sum()), all_on=int((self.n_off[m] == 0).sum()),
                    expected_all_off=float((w * p_off).sum()), expected_all_on=float((w * p_on).sum()))
