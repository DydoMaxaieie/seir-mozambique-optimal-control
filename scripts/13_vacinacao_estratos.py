"""Extensão: doses, cobertura e distanciamento médio por região e grupo etário (estratégia combinada, caso base)."""
import sys, pathlib, pickle, csv
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from seirmoz.model import N
from seirmoz.paths import RESULTS
x0, beta, kappa, out, J0 = pickle.load(open(RESULTS / "controlo" / "base.pkl", "rb"))[:5]
v, u, X, J = out["vu"]; Xs = out["sem"][2]
REG = ["Norte", "Centro", "Sul-interior", "Sul-metropolitano"]; GR = ["0-4", "5-17", "18-39", "40-59", "60+"]
doses = X[-1][12] - x0[12]; cob = doses / N; um = u.mean(0)
d0 = x0[10] + x0[11]; ob_c = X[-1][10] + X[-1][11] - d0; ob_s = Xs[-1][10] + Xs[-1][11] - d0
np.set_printoptions(suppress=True, linewidth=150)
print("doses (mil)\n", np.round(doses / 1e3, 1)); print("cobertura %\n", np.round(100 * cob, 1)); print("u medio\n", np.round(um, 3))
print("obitos sem\n", np.round(ob_s)); print("obitos com\n", np.round(ob_c))
print("total doses", doses.sum(), "obitos por idade sem", np.round(ob_s.sum(0)), "com", np.round(ob_c.sum(0)))
with open(RESULTS / "controlo" / "estratos.csv", "w", newline="", encoding="utf8") as f:
    w = csv.writer(f); w.writerow(["regiao", "grupo", "populacao", "doses", "cobertura_pct", "u_medio", "obitos_sem", "obitos_com"])
    for k in range(4):
        for g in range(5):
            w.writerow([REG[k], GR[g], round(N[k, g]), round(doses[k, g]), round(100 * cob[k, g], 2), round(um[k, g], 3), round(ob_s[k, g], 1), round(ob_c[k, g], 1)])
