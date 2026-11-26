"""Extensão: validação fora da amostra por origem móvel.

Em cada origem D ajusta-se beta só com os dias < D (troços já iniciados) e prevê-se H dias à frente
com o último beta ajustado. Compara-se o RMSE da previsão com o da persistência (média móvel dos últimos 7 dias).
Os pontos de mudança de beta são os da dissertação (ficam fixos, estimados com a série completa).
"""
import sys, pathlib, csv
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.optimize import least_squares
from seirmoz import calibration as cal
from seirmoz.paths import RESULTS

H = 14
ORIGENS = ["2021-03-15", "2021-04-15", "2021-05-15", "2021-06-10"]
Tfull = cal.T
b_full, a_full = cal.load_params(0.1)
rho = 0.1
rows = []
for D_s in ORIGENS:
    D = cal.dates.index(D_s)
    cal.T = D + H
    nseg = [1 + sum(cal.dates.index(c) < D for c in cal.cuts[r]) for r in cal.regs]
    b = b_full.copy(); a = a_full.copy()
    for it in range(2):
        for k in range(cal.K):
            m = nseg[k]
            def res(p):
                bb = b.copy(); bb[k, :m] = np.exp(p[:m]); bb[k, m:] = np.exp(p[m - 1])
                aa = a.copy(); aa[k] = np.exp(p[m])
                z = cal.sim(bb, rho, a=aa)
                return np.sqrt(z[:D, k]) - np.sqrt(cal.ys[:D, k])
            starts = [np.r_[np.log(b[k, :m]), np.log(a[k])]] + [np.r_[np.log(np.full(m, s)), 0.0] for s in (0.004, 0.01, 0.025)]
            lo = np.r_[np.log(np.full(m, 5e-4)), np.log(0.05)]; hi = np.r_[np.log(np.full(m, 0.5)), np.log(20)]
            best = min((least_squares(res, p0, bounds=(lo, hi)) for p0 in starts), key=lambda r: r.cost)
            b[k, :m] = np.exp(best.x[:m]); b[k, m:] = np.exp(best.x[m - 1]); a[k] = np.exp(best.x[m])
    z = cal.sim(b, rho, a=a)
    obs = cal.ys[D:D + H]; prev = z[D:D + H]
    pers = np.tile(cal.ys[D - 7:D].mean(0), (H, 1))
    ins = np.sqrt(((z[:D] - cal.ys[:D]) ** 2).mean(0))
    rm = np.sqrt(((prev - obs) ** 2).mean(0)); rp = np.sqrt(((pers - obs) ** 2).mean(0))
    # erro relativo ao nível médio observado no horizonte
    nivel = obs.mean(0)
    for k, r in enumerate(cal.regs):
        rows.append([D_s, r, round(ins[k], 1), round(rm[k], 1), round(rp[k], 1), round(float(nivel[k]), 1)])
    print(D_s, "modelo", np.round(rm, 1), "persistência", np.round(rp, 1), "nível", np.round(nivel, 1), flush=True)
cal.T = Tfull
with open(RESULTS / "controlo" / "validacao_previsao.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["origem", "macrorregiao", "rmse_ajuste", "rmse_previsao_14d", "rmse_persistencia_14d", "nivel_medio_observado"]); w.writerows(rows)
