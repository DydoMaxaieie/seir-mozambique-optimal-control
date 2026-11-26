"""Extensão: 40 perturbações admissíveis do controlo óptimo com a normalização por idade (J deve aumentar)."""
import sys, pathlib, pickle
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from seirmoz import calibration as cal, model
cal.set_norm('idade'); model.set_norm('idade')
from seirmoz.control import setup, K, G
P, J0, k, b = setup(0.1)
v, u, X, J, h, it = P.solve()
print("J óptimo", J, "iterações", it)
rng = np.random.default_rng(1); inc = []
for i in range(40):
    dv = rng.standard_normal(v.shape) * P.vmax.max() * 0.02; du = rng.standard_normal(u.shape) * 0.02
    v2 = np.clip(v + dv, 0, P.vmax[None]); u2 = np.clip(u + du, 0, P.umax)
    inc.append(P.forward(v2, u2)[1] - J)
print("perturbações: mínimo ΔJ", min(inc), " nº negativos", sum(x < 0 for x in inc))
