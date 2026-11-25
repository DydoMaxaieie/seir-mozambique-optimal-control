"""Capítulo 5: calibração por mínimos quadrados (por rho), métricas e R0 por bloco.

Uso: python scripts/02_calibracao.py [--refazer]
Sem --refazer, usa os parâmetros guardados em results/calibracao (os da dissertação).
Com --refazer, recalcula de raiz com rho=0.1 (cerca de 2 min) e compara.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from seirmoz import calibration as cal

rho = 0.1
if "--refazer" in sys.argv:
    beta, a = cal.fit(rho, passes=3)[:2]
    ref = cal.load_params(rho)[0]
    print("diferença máxima de beta face aos guardados:", np.abs(beta - ref).max())
else:
    beta, a = cal.load_params(rho)[:2]
np.set_printoptions(precision=4, suppress=True)
print("beta por região e bloco:\n", beta)
z = cal.sim(beta, rho, a)
print("RMSE, R2:", cal.metrics(z))
print("R0 por bloco:")
for k in range(4):
    print(" ", [round(float(cal.R0_block(k, beta[k, b])), 2) for b in range(beta.shape[1])])
