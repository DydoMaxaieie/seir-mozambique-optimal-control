"""Extensão: quatro estratégias com a convenção padrão C_gh I_h / N_h (comparação com o caso base)."""
import sys, pathlib, pickle, csv
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from seirmoz import calibration as cal, model
mode = sys.argv[1] if len(sys.argv) > 1 else 'idade'
cal.set_norm(mode); model.set_norm(mode)
import importlib.util
spec = importlib.util.spec_from_file_location("c3", pathlib.Path(__file__).parent / "03_controlo.py"); c3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c3)
from seirmoz.control import K, G, N
from seirmoz.paths import RESULTS
OUT = RESULTS / ("controlo" if mode == 'total' else "controlo_idade"); OUT.mkdir(exist_ok=True)
P, J0, k, out = c3.run(0.1, 0.5, 0.1)
print("kappa", round(k, 3), "J0", round(J0, 1), flush=True)
res = {}
for s, nome in (("sem", "Sem controlo"), ("v", "Só vacinação"), ("u", "Só distanciamento"), ("vu", "Vacinação e distanciamento")):
    v, u, X, J = out[s]; m = c3.metrics(P, X, v, u, J)
    ub = (u.mean(0) * N).sum(0) / N.sum(0) if s in ("u", "vu") else np.zeros(G)
    ur = (u.mean(0) * N).sum(1) / N.sum(1) if s in ("u", "vu") else np.zeros(K)
    cov = X[-1][12] / N
    covg = (X[-1][12].sum(0)) / N.sum(0)
    res[s] = dict(m, u_idade=ub, u_regiao=ur, cob_idade=covg)
    print(f"{nome:28s} óbitos {m['obitos']:9.1f} infecções {m['infeccoes']/1e6:5.2f}M doses {m['doses']/1e3:7.0f}k J/J0 {J/J0:.3f}")
    print("   óbitos região", np.round(m['obitos_regiao']), "idade", np.round(m['obitos_idade']))
    print("   u por região", np.round(ur, 3), "u por idade", np.round(ub, 3))
    print("   cobertura por idade", np.round(covg, 4), flush=True)
pickle.dump((P.x0, P.beta, P.kappa, out, J0, res), open(OUT / "base.pkl", "wb"))
