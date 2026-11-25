"""Capítulo 3: verificação do FBS no problema escalar (Tabela 3.1) e no SEIR homogéneo (Tabela 3.2)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
from seirmoz import verif_escalar as lq, verif_seir as seir

print("== Problema linear-quadrático escalar (Tabela 3.1) ==")
for m in ("euler", "rk4avg", "hermite"):
    errs = []
    for n in (10, 20, 40, 80, 160):
        t, x, q, u, it = lq.fbs(m, n)
        errs.append(np.abs(x - lq.xe(t)).max())
    ordem = [np.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
    print(m, ["%.2e" % e for e in errs], np.round(ordem, 2))
t, x, q, u, it = lq.fbs("hermite", 80)
print("q(0) numérico", q[0], "exacto", lq.qe(0))

print("== SEIR homogéneo (Tabela 3.2) ==")
seir.main()
