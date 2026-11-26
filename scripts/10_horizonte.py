"""Extensão: redução de óbitos em função do horizonte, com e sem custo terminal, e óbitos a 600 dias."""
import sys, pathlib, csv
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from seirmoz.extensoes import *
from seirmoz.model import step, rhs
from seirmoz.paths import RESULTS

TEND = 600
def tail(P, x, T):
    """Continua a integração sem controlo até TEND e devolve os óbitos acumulados."""
    n = int(round((TEND - T) / P.h)); z = np.zeros((K, G))
    for _ in range(n):
        x = step(rhs, x, P.h, z, z, P.beta, P.kappa)
    return (x[10] + x[11]).sum() - (P.x0[10] + P.x0[11]).sum()

rows = []
for T in (120, 180, 240, 300):
    for term in (False, True):
        P, J0, k, b = setup(0.1, T=T); P.__class__ = ProbT; P.set_terminal(term)
        z = np.zeros((P.n, K, G)); X0, _ = P.forward(z, z)
        d0 = P.deaths(X0)[0]; d0_end = tail(P, X0[-1], T)
        v, u, X, J, h, it = P.solve(it=1500, relax=0.5 if T <= 180 else 0.25)
        dT = P.deaths(X)[0]; dend = tail(P, X[-1], T)
        rows.append([T, int(term), round(d0), round(dT), round(100 * (1 - dT / d0), 1), round(d0_end), round(dend), round(100 * (1 - dend / d0_end), 1), it])
        print(rows[-1], flush=True)
with open(RESULTS / "controlo" / "horizonte.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["T", "custo_terminal", "obitos_sem_T", "obitos_vu_T", "reducao_T_pct", "obitos_sem_600", "obitos_vu_600", "reducao_600_pct", "iteracoes"]); w.writerows(rows)
