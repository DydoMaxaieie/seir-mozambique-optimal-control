"""Capítulo 6: quatro estratégias (Tabelas 6.2 e 6.3) e sensibilidades (Tabelas 6.4 e 6.5).

Uso: python scripts/03_controlo.py [--rapido]   (--rapido só corre o caso base)
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import csv, pickle, sys
import numpy as np
from seirmoz.control import setup, K, G, N, TE
from seirmoz.paths import RESULTS

OUT = RESULTS / "controlo"
OUT.mkdir(parents=True, exist_ok=True)


def metrics(P, X, v, u, J):
    d0 = P.x0[10] + P.x0[11]; dT = X[-1][10] + X[-1][11]
    inc = np.array([X[i][3] / TE for i in range(len(X))])  # fluxo E -> I
    return dict(obitos=(dT - d0).sum(), infeccoes=np.trapezoid(inc.sum((1, 2)), dx=P.h),
                doses=X[-1][12].sum(), J=J, cobertura_max=(X[-1][12] / N).max(),
                obitos_regiao=(dT - d0).sum(1), obitos_idade=(dT - d0).sum(0))


def run(rho, eta_u, eta_v, T=120):
    P, J0, k, b = setup(rho, T=T, eta_u=eta_u, eta_v=eta_v)
    z = np.zeros((P.n, K, G)); out = {}
    X0, Jz = P.forward(z, z); out["sem"] = (z, z, X0, Jz)
    v, u, X, J, h, it = P.solve(free_u=False); out["v"] = (v, u, X, J)
    v, u, X, J, h, it = P.solve(free_v=False); out["u"] = (v, u, X, J)
    v, u, X, J, h, it = P.solve(); out["vu"] = (v, u, X, J)
    return P, J0, k, out


if __name__ == "__main__":
    P, J0, k, out = run(0.1, 0.5, 0.1)
    pickle.dump((P.x0, P.beta, P.kappa, out, J0), open(OUT / "base.pkl", "wb"))
    print("kappa", round(k, 3), "J0", round(J0, 1))
    rows = []
    for s, nome in (("sem", "Sem controlo"), ("v", "Só vacinação"), ("u", "Só distanciamento"), ("vu", "Vacinação e distanciamento")):
        v, u, X, J = out[s]; m = metrics(P, X, v, u, J)
        print(f"{nome:28s} óbitos {m['obitos']:9.1f}  infecções {m['infeccoes']/1e6:5.2f} M  doses {m['doses']/1e3:7.0f} mil  J/J0 {J/J0:.3f}")
        rows.append([nome, round(m["obitos"], 1), round(m["infeccoes"]), round(m["doses"]), round(J / J0, 4)] + [round(x) for x in m["obitos_regiao"]])
    with open(OUT / "estrategias_base.csv", "w", newline="", encoding="utf8") as f:
        w = csv.writer(f); w.writerow(["estrategia", "obitos", "infeccoes", "doses", "J_sobre_J0", "Norte", "Centro", "Sul-interior", "Sul-metropolitano"]); w.writerows(rows)
    if "--rapido" in sys.argv:
        sys.exit()
    rows = []
    for rho in (0.05, 0.1, 0.2):
        P, J0, k, out = run(rho, 0.5, 0.1)
        for s in ("sem", "vu"):
            v, u, X, J = out[s]; m = metrics(P, X, v, u, J)
            rows.append(("rho", rho, 0.5, 0.1, s, round(m["obitos"], 1), round(m["infeccoes"]), round(m["doses"]), round(J / J0, 4), round(float(u.mean()), 4), round(k, 3)))
    for eu in (0.1, 0.25, 0.5, 1.0):
        for ev in (0.05, 0.1, 0.5):
            P, J0, k, out = run(0.1, eu, ev)
            v, u, X, J = out["vu"]; m = metrics(P, X, v, u, J)
            rows.append(("pesos", 0.1, eu, ev, "vu", round(m["obitos"], 1), round(m["infeccoes"]), round(m["doses"]), round(J / J0, 4), round(float(u.mean()), 4), round(k, 3)))
    with open(OUT / "sensibilidade.csv", "w", newline="", encoding="utf8") as f:
        w = csv.writer(f); w.writerow(["tipo", "rho", "eta_u", "eta_v", "cenario", "obitos", "infeccoes", "doses", "J_sobre_J0", "u_medio", "kappa"]); w.writerows(rows)
    print("escrito", OUT / "sensibilidade.csv")
