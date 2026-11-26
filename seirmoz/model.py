"""Dinâmica de 13 compartimentos por estrato (4.3)-(4.15) e sistema adjunto (4.43)-(4.53)."""
import numpy as np
from . import inputs

_d = inputs.load()
N, C, Th, ag = _d["N"], _d["C"], _d["Th"], _d["ag"]
K, G = N.shape
TE, TI, TQ, THW, THC, TV, EFF, OMEGA = 5.2, 10.0, 7.0, 8.0, 12.0, 14.0, 0.70, 0.5
Nhat = Th.T @ N.sum(1)
Nhat_g = Th.T @ N  # população presente por região e grupo etário (K x G)
NORM = 'total'  # 'total': denominador N̂_m (dissertação); 'idade': N̂_{mh} (convenção padrão)


def set_norm(mode):
    global NORM
    assert mode in ('total', 'idade')
    NORM = mode
W = (Th / Nhat[None, :]) @ Th.T  # W_kl = sum_m th_km th_lm / Nhat_m
IDX = dict(Su=0, Sv=1, Sp=2, E=3, I=4, Q=5, Hw=6, Hc=7, R=8, RH=9, D=10, DH=11, Vw=12)

# ---- parâmetros clínicos por grupo (Tabelas 6.1 e D.1) ----
# intervalos 0-9, 10-19, ..., 80+ (Verity et al., 2020: IHR; Ferguson et al., 2020: p^c e IFR)
IHR = np.array([0, 0.0408, 1.04, 3.43, 4.25, 8.16, 11.8, 16.6, 18.4]) / 100
PC = np.array([5, 5, 5, 5, 6.3, 12.2, 27.4, 43.2, 70.9]) / 100
IFR = np.array([0.002, 0.006, 0.03, 0.08, 0.15, 0.6, 2.2, 5.1, 9.3]) / 100
Wt = np.zeros((G, 9))
Wt[0, 0] = 1; Wt[1, 0] = 5; Wt[1, 1] = 8; Wt[2, 1] = 2; Wt[2, 2] = 10; Wt[2, 3] = 10
Wt[3, 4] = 10; Wt[3, 5] = 10
Wt[4, 6] = .60; Wt[4, 7] = .29; Wt[4, 8] = .11
Wt = Wt / Wt.sum(1, keepdims=True)
ph = np.maximum(Wt @ IHR, 1e-4)
ifr = Wt @ IFR
pc = (Wt @ (IHR * PC)) / np.maximum(Wt @ IHR, 1e-12)
pc[0] = 0.05
PHI_OUT = 0.30  # fracção dos óbitos fora do hospital
fq = PHI_OUT * ifr / (1 - ph)
h_ = (1 - PHI_OUT) * ifr / ph  # probabilidade de morrer dado internado
fw = 0.3 * h_
fc = np.minimum((h_ - fw) / ((1 - fw) * pc), 0.5)
mu_q = fq / ((1 - fq) * TQ)
mu_w = fw / ((1 - fw) * THW)
mu_c = fc / ((1 - fc) * THC)


def lam_f(I, u, beta):
    prev = (Th.T @ I) / (Nhat_g if NORM == 'idade' else Nhat[:, None])
    return (1 - u) * beta[:, None] * (Th @ (prev @ C.T))


def rhs(x, v, u, beta, kappa):
    Su, Sv, Sp, E, I, Q, Hw, Hc, R, RH, Dd, DH, Vw = x
    lam = lam_f(I, u, beta)
    dx = np.empty_like(x)
    dx[0] = -(lam + v) * Su
    dx[1] = v * Su - (lam + 1 / TV) * Sv
    dx[2] = (1 - EFF) / TV * Sv - (1 - OMEGA) * lam * Sp
    dx[3] = lam * (Su + Sv + (1 - OMEGA) * Sp) - E / TE
    dx[4] = E / TE - I / TI
    dx[5] = (1 - ph)[None, :] / TI * I - (1 / TQ + kappa * mu_q)[None, :] * Q
    dx[6] = ph[None, :] / TI * I - (1 / THW + kappa * mu_w)[None, :] * Hw
    dx[7] = pc[None, :] / THW * Hw - (1 / THC + kappa * mu_c)[None, :] * Hc
    dx[8] = EFF / TV * Sv + Q / TQ
    dx[9] = (1 - pc)[None, :] / THW * Hw + Hc / THC
    dx[10] = kappa * mu_q[None, :] * Q
    dx[11] = kappa * (mu_w[None, :] * Hw + mu_c[None, :] * Hc)
    dx[12] = v * Su
    return dx


def adj_rhs(L, x, v, u, beta, kappa, A, AH):
    Su, Sv, Sp, E, I, Q, Hw, Hc, R, RH, Dd, DH, Vw = x
    LSu, LSv, LSp, LE, LI, LQ, LHw, LHc, LR, LRH, LD, LDH, LV = L
    lam = lam_f(I, u, beta)
    Psi = (LE - LSu) * Su + (LE - LSv) * Sv + (1 - OMEGA) * (LE - LSp) * Sp
    M = Psi * (1 - u) * beta[:, None]
    if NORM == 'idade':
        Xi = Th @ ((Th.T @ (M @ C)) / Nhat_g)
    else:
        Xi = W.T @ (M @ C)  # Xi_lh, eq. (4.48)
    dL = np.empty_like(L)
    dL[0] = (lam + v) * LSu - lam * LE - v * LSv - v * LV
    dL[1] = (lam + 1 / TV) * LSv - lam * LE - (1 - EFF) / TV * LSp - EFF / TV * LR
    dL[2] = (1 - OMEGA) * lam * (LSp - LE)
    dL[3] = (LE - LI) / TE
    dL[4] = (LI - (1 - ph)[None, :] * LQ - ph[None, :] * LHw) / TI - Xi
    dL[5] = (1 / TQ + kappa * mu_q)[None, :] * LQ - LR / TQ - kappa * mu_q[None, :] * LD
    dL[6] = (1 / THW + kappa * mu_w)[None, :] * LHw - pc[None, :] / THW * LHc - (1 - pc)[None, :] / THW * LRH - kappa * mu_w[None, :] * LDH
    dL[7] = (1 / THC + kappa * mu_c)[None, :] * LHc - LRH / THC - kappa * mu_c[None, :] * LDH
    dL[8] = 0; dL[9] = 0; dL[12] = 0
    dL[10] = -A; dL[11] = -AH
    return dL, Psi, lam


def step(f, x, h, *a):
    """Passo de Runge-Kutta de quarta ordem."""
    k1 = f(x, *a); k2 = f(x + h / 2 * k1, *a); k3 = f(x + h / 2 * k2, *a); k4 = f(x + h * k3, *a)
    return x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def window_state(rho, kappa, h=0.25):
    """Integra a janela de calibração (u = v = 0) e devolve o estado em 30/06/2021."""
    from . import calibration as cal
    b, a = cal.load_params(rho)
    z0 = cal.ys[:7].mean(0) * a
    I0 = (z0 * TI / rho)[:, None] * ag[None, :]
    E0 = (z0 * TE / rho)[:, None] * ag[None, :]
    x = np.zeros((13, K, G)); x[0] = N - I0 - E0; x[3] = E0; x[4] = I0
    n = int(round(182 / h)); zs = []
    v = np.zeros((K, G)); u = np.zeros((K, G))
    for i in range(n):
        t = i * h
        bt = b[np.arange(K), cal.seg[:, min(int(t), 181)]]
        if i % int(1 / h) == 0:
            zs.append(rho * x[3].sum(1) / TE)
        x = step(rhs, x, h, v, u, bt, kappa)
    return x, np.array(zs), b, cal.ys
