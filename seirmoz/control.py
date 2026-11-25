"""Problema de controlo óptimo e varrimento progressivo-regressivo (Capítulos 4 e 6)."""
import numpy as np
from scipy.optimize import brentq
from .model import *
from .model import window_state, adj_rhs, rhs, step

OBITOS_JANELA = 712  # óbitos notificados entre 31/12/2020 (166) e 30/06/2021 (878)


def kappa_rho(rho):
    """Modulador de mortalidade que reproduz os 712 óbitos da janela."""
    f = lambda k: (lambda x: (x[10] + x[11]).sum() - OBITOS_JANELA)(window_state(rho, k)[0])
    return brentq(f, 0.3, 3.0, xtol=1e-4)


class Prob:
    def __init__(s, x0, beta, kappa, T=180, h=0.5, A=1.0, AH=1.0, B=1.0, Cw=1.0, vmax=None, umax=0.5, vbar=0.6):
        s.x0 = x0; s.beta = beta; s.kappa = kappa; s.T = T; s.h = h; s.n = int(round(T / h))
        s.A = np.full((K, G), A); s.AH = np.full((K, G), AH); s.B = np.full((K, G), B); s.Cw = np.full((K, G), Cw)
        s.vmax = np.array([0, 0, .005, .005, .005])[None, :].repeat(K, 0) if vmax is None else vmax
        s.umax = umax; s.vbar = vbar

    def forward(s, v, u):
        x = s.x0.copy(); X = [x]; J = 0.0
        def cost(x, v, u):
            return (s.A * x[10] + s.AH * x[11] + 0.5 * s.B * v ** 2 + 0.5 * s.Cw * N * u ** 2).sum()
        for i in range(s.n):
            vi = np.where(x[12] >= s.vbar * N, 0.0, v[i])  # restrição de cobertura
            c0 = cost(x, vi, u[i])
            x = step(rhs, x, s.h, vi, u[i], s.beta, s.kappa); X.append(x)
            J += 0.5 * s.h * (c0 + cost(x, vi, u[i]))
        return np.array(X), J

    def backward(s, X, v, u):
        L = np.zeros((13, K, G)); Ls = [L]
        for i in range(s.n - 1, -1, -1):
            xa, xb = X[i], X[i + 1]; xm = (xa + xb) / 2
            f = lambda L_, x_: adj_rhs(L_, x_, v[i], u[i], s.beta, s.kappa, s.A, s.AH)[0]
            h = -s.h
            k1 = f(L, xb); k2 = f(L + h / 2 * k1, xm); k3 = f(L + h / 2 * k2, xm); k4 = f(L + h * k3, xa)
            L = L + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4); Ls.append(L)
        return np.array(Ls[::-1])

    def characterise(s, X, Ls, v, u):
        """Controlos (4.54) e (4.55) a partir do estado e dos co-estados."""
        vn = np.zeros_like(v); un = np.zeros_like(u)
        for i in range(s.n):
            x = (X[i] + X[i + 1]) / 2; L = (Ls[i] + Ls[i + 1]) / 2
            _, Psi, lam = adj_rhs(L, x, v[i], u[i], s.beta, s.kappa, s.A, s.AH)
            lt = lam / (1 - u[i])
            vn[i] = np.clip((L[0] - L[1] - L[12]) * x[0] / s.B, 0, s.vmax)
            un[i] = np.clip(lt * Psi / (s.Cw * N), 0, s.umax)
        return vn, un

    def solve(s, v0=None, u0=None, it=300, relax=0.5, tol=1e-5, verbose=False, free_v=True, free_u=True):
        v = np.zeros((s.n, K, G)) if v0 is None else v0.copy()
        u = np.zeros((s.n, K, G)) if u0 is None else u0.copy()
        hist = []
        for k in range(it):
            X, J = s.forward(v, u); Ls = s.backward(X, v, u); vn, un = s.characterise(X, Ls, v, u)
            if not free_v: vn[:] = 0
            if not free_u: un[:] = 0
            dv = np.abs(vn - v).max() / max(s.vmax.max(), 1e-12); du = np.abs(un - u).max() / s.umax
            v = relax * vn + (1 - relax) * v; u = relax * un + (1 - relax) * u; hist.append(J)
            if verbose and k % 10 == 0:
                print(k, round(J, 1), round(dv, 5), round(du, 5))
            if max(dv, du) < tol:
                break
        X, J = s.forward(v, u)
        return v, u, X, J, hist, k + 1


def setup(rho, T=120, eta_u=0.5, eta_v=0.1, umax=0.5, vmaxv=0.005, vbar=0.6, kap=None):
    """Cenário base: pesos B = b N_kg e C = c calibrados por eta_v e eta_u (eq. 6.3)."""
    k = kappa_rho(rho) if kap is None else kap
    x, z, b, ys = window_state(rho, k)
    P = Prob(x, b[:, 3], k, T=T, umax=umax, vbar=vbar)
    P.vmax = np.array([0, 0, 1, 1, 1])[None, :].repeat(K, 0) * vmaxv
    z0 = np.zeros((P.n, K, G)); X, J0 = P.forward(z0, z0)
    elig = (P.vmax > 0)
    b0 = eta_v * J0 / (0.5 * (N * elig * P.vmax ** 2).sum() * T)
    P.B = b0 * N
    c0 = eta_u * J0 / (0.5 * N.sum() * umax ** 2 * T); P.Cw = np.full((K, G), c0)
    return P, J0, k, b
