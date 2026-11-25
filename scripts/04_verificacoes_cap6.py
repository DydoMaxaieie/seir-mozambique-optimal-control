"""Capítulo 6: gradiente das adjuntas vs diferenças finitas, perturbações do óptimo e convergência em h."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
from seirmoz.control import setup, adj_rhs, K, G, N

P, J0, k, b = setup(0.1)
rng = np.random.default_rng(1)
v = np.zeros((P.n, K, G)); u = np.full((P.n, K, G), 0.15); v[:] = 0.002 * (P.vmax > 0)
X, J = P.forward(v, u); Ls = P.backward(X, v, u)
gu = np.zeros_like(u); gv = np.zeros_like(v)
for i in range(P.n):
    x = (X[i] + X[i + 1]) / 2; L = (Ls[i] + Ls[i + 1]) / 2
    _, Psi, lam = adj_rhs(L, x, v[i], u[i], P.beta, P.kappa, P.A, P.AH)
    gu[i] = P.h * (P.Cw * N * u[i] - lam / (1 - u[i]) * Psi)
    gv[i] = P.h * (P.B * v[i] - (L[0] - L[1] - L[12]) * x[0])
print("1) gradiente (adjunta vs diferenças finitas centradas)")
erros = []
for (i, kk, g, tipo) in [(10, 1, 2, "u"), (40, 3, 3, "u"), (70, 0, 4, "u"), (100, 2, 2, "u"), (20, 1, 3, "v"), (60, 3, 4, "v"), (90, 0, 4, "v")]:
    eps = 1e-4 if tipo == "u" else 1e-6
    vp, up, vm, um = v.copy(), u.copy(), v.copy(), u.copy()
    if tipo == "u": up[i, kk, g] += eps; um[i, kk, g] -= eps
    else: vp[i, kk, g] += eps; vm[i, kk, g] -= eps
    fd = (P.forward(vp, up)[1] - P.forward(vm, um)[1]) / (2 * eps)
    an = (gu if tipo == "u" else gv)[i, kk, g]
    erros.append(abs(fd - an) / max(abs(fd), 1e-9))
    print(f"  {tipo} (t={i}, k={kk}, g={g}) FD {fd:.3f} adjunta {an:.3f} erro rel {erros[-1]:.1e}")
print("  erro relativo máximo", f"{max(erros):.1e}")

print("2) perturbações aleatórias admissíveis do óptimo")
v, u, X, J, h, it = P.solve(it=200)
piores = 0
for r in range(40):
    dv = rng.normal(size=v.shape) * P.vmax.max() * 0.1 * (P.vmax > 0); du = rng.normal(size=u.shape) * 0.05
    Jp = P.forward(np.clip(v + dv, 0, P.vmax), np.clip(u + du, 0, P.umax))[1]
    piores += Jp >= J
print(f"  {piores}/40 perturbações aumentaram J; iterações do varrimento: {it}")

print("3) convergência em h")
for hh in (1.0, 0.5, 0.25):
    P, J0, k, b = setup(0.1)
    P.h = hh; P.n = int(round(P.T / hh))
    v, u, X, J, _, it = P.solve()
    d = (X[-1][10] + X[-1][11]).sum() - (P.x0[10] + P.x0[11]).sum()
    print(f"  h = {hh}: óbitos {d:.1f}  iterações {it}")
