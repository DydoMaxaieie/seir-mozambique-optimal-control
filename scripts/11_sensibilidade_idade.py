"""Extensão: sensibilidade a rho e aos pesos com a normalização por idade."""
import sys, pathlib, csv, importlib.util
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from seirmoz import calibration as cal, model
cal.set_norm('idade'); model.set_norm('idade')
spec = importlib.util.spec_from_file_location("c3", pathlib.Path(__file__).parent / "03_controlo.py"); c3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c3)
from seirmoz.paths import RESULTS
OUT = RESULTS / "controlo_idade"; OUT.mkdir(exist_ok=True)
rows = []
for rho in (0.05, 0.1, 0.2):
    P, J0, k, out = c3.run(rho, 0.5, 0.1)
    for s in ("sem", "v", "u", "vu"):
        v, u, X, J = out[s]; m = c3.metrics(P, X, v, u, J)
        rows.append(("rho", rho, 0.5, 0.1, s, round(m["obitos"], 1), round(m["infeccoes"]), round(m["doses"]), round(J / J0, 4), round(float(u.mean()), 4), round(k, 3)))
    print(rho, [r[5] for r in rows[-4:]], flush=True)
for eu in (0.1, 0.25, 0.5, 1.0):
    for ev in (0.05, 0.1, 0.5):
        P, J0, k, out = c3.run(0.1, eu, ev)
        v, u, X, J = out["vu"]; m = c3.metrics(P, X, v, u, J)
        rows.append(("pesos", 0.1, eu, ev, "vu", round(m["obitos"], 1), round(m["infeccoes"]), round(m["doses"]), round(J / J0, 4), round(float(u.mean()), 4), round(k, 3)))
        print(eu, ev, round(m["obitos"]), flush=True)
with open(OUT / "sensibilidade.csv", "w", newline="", encoding="utf8") as f:
    w = csv.writer(f); w.writerow(["tipo", "rho", "eta_u", "eta_v", "cenario", "obitos", "infeccoes", "doses", "J_sobre_J0", "u_medio", "kappa"]); w.writerows(rows)
