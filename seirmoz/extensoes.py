"""Extensão: custo terminal pelos óbitos comprometidos (compartimentos E, I, Q, H^w, H^c) no fim do horizonte."""
import numpy as np
from . import model
from .model import K, G, N, IDX, TI, TQ, THW, THC, ph, pc, mu_q, mu_w, mu_c, adj_rhs, rhs, step
from .control import Prob, setup


def committed(kappa):
    """Probabilidade de morte futura, sem nova transmissão, a partir de cada compartimento (vector por grupo)."""
    dQ = kappa * mu_q / (1 / TQ + kappa * mu_q)
    dc = kappa * mu_c / (1 / THC + kappa * mu_c)
    dw = (kappa * mu_w + pc / THW * dc) / (1 / THW + kappa * mu_w)
    dI = (1 - ph) * dQ + ph * dw
    return dQ, dw, dc, dI


class ProbT(Prob):
    def set_terminal(s, on=True):
        dQ, dw, dc, dI = committed(s.kappa)
        t = np.zeros((13, K, G))
        if on:
            t[IDX['E']] = s.AH[None, 0] * dI if False else (dI[None, :] * s.A)  # A = AH = 1 no caso base
            t[IDX['I']] = dI[None, :] * s.A
            t[IDX['Q']] = dQ[None, :] * s.A
            t[IDX['Hw']] = dw[None, :] * s.AH
            t[IDX['Hc']] = dc[None, :] * s.AH
        s.term = t

    def forward(s, v, u):
        X, J = Prob.forward(s, v, u)
        return X, J + float((s.term * X[-1]).sum())

    def backward(s, X, v, u):
        L = s.term.copy(); Ls = [L]
        for i in range(s.n - 1, -1, -1):
            xa, xb = X[i], X[i + 1]; xm = (xa + xb) / 2
            f = lambda L_, x_: adj_rhs(L_, x_, v[i], u[i], s.beta, s.kappa, s.A, s.AH)[0]
            h = -s.h
            k1 = f(L, xb); k2 = f(L + h / 2 * k1, xm); k3 = f(L + h / 2 * k2, xm); k4 = f(L + h * k3, xa)
            L = L + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4); Ls.append(L)
        return np.array(Ls[::-1])

    def deaths(s, X):
        d = (X[-1][10] + X[-1][11]).sum() - (s.x0[10] + s.x0[11]).sum()
        c = float((s.term * X[-1]).sum()) if hasattr(s, 'term') else 0.0
        return d, c
