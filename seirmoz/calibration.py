"""Calibração de beta_k^(j) por macrorregião e troço (Secção 5.2).

O sistema reduzido (S, E, I) integra-se com LSODA; ajusta-se a raiz quadrada da
incidência esperada à média móvel de sete dias, uma macrorregião de cada vez.
"""
import csv
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from . import inputs
from .paths import DATA, CALIB

_d = inputs.load()
N, C, Th, ag = _d["N"], _d["C"], _d["Th"], _d["ag"]
K, G = N.shape
TE, TI = 5.2, 10.0
regs = inputs.REGS

# ---- casos ----
with open(DATA / "casos_macrorregiao.csv", encoding="utf8") as _f:
    _R = list(csv.DictReader(_f))
dates = [r["data"] for r in _R]
y = np.array([[float(r[k]) for k in regs] for r in _R])
T = len(y)


def ma7(x):
    p = np.pad(x, (3, 3), mode="edge")
    return np.convolve(p, np.ones(7) / 7, "valid")


ys = np.column_stack([ma7(y[:, k]) for k in range(K)])

# datas de mudança de beta (Tabela 5.5), obtidas por regressão linear por troços (programação dinâmica)
cuts = {
    "Norte": ["2021-02-26", "2021-04-22", "2021-06-05"],
    "Centro": ["2021-02-08", "2021-04-17", "2021-05-26"],
    "Sul-interior": ["2021-01-21", "2021-03-05", "2021-05-08"],
    "Sul-metropolitano": ["2021-01-23", "2021-05-19", "2021-06-09"],
}
idx = lambda ds: [dates.index(x) for x in ds]
seg = np.zeros((K, T), int)
for k, r in enumerate(regs):
    for c in idx(cuts[r]):
        seg[k, c:] += 1
Nhat = Th.T @ N.sum(1)  # população presente em cada região


def sim(beta, rho, a=None, imm0=None):
    """Casos confirmados esperados por dia (T x K). beta: K x 4; a: factor de correcção do estado inicial."""
    z0 = ys[:7].mean(0) * (np.ones(K) if a is None else a)
    I0 = (z0 * TI / rho)[:, None] * ag[None, :]
    E0 = (z0 * TE / rho)[:, None] * ag[None, :]
    S0 = N - I0 - E0
    if imm0 is not None:
        S0 = S0 - imm0
    x0 = np.r_[S0.ravel(), E0.ravel(), I0.ravel()]
    n = K * G

    def f(t, x):
        S = x[:n].reshape(K, G); E = x[n:2 * n].reshape(K, G); I = x[2 * n:].reshape(K, G)
        b = beta[np.arange(K), seg[:, min(int(t), T - 1)]]
        prev = (Th.T @ I) / Nhat[:, None]
        lam = b[:, None] * (Th @ (prev @ C.T))
        new = lam * S
        return np.r_[(-new).ravel(), (new - E / TE).ravel(), (E / TE - I / TI).ravel()]

    sol = solve_ivp(f, (0, T), x0, t_eval=np.arange(T) + 0.5, method="LSODA", rtol=1e-7, atol=1e-3)
    E = sol.y[n:2 * n].reshape(K, G, -1).sum(1)
    return (rho * E / TE).T


cost = lambda z: 0.5 * ((np.sqrt(z) - np.sqrt(ys)) ** 2).sum()


def fit(rho, b=None, a=None, passes=3, verbose=True):
    """Ajuste por coordenadas: uma macrorregião de cada vez, com quatro pontos de partida."""
    b = np.full((K, 4), 0.012) if b is None else b.copy()
    a = np.ones(K) if a is None else a.copy()
    for it in range(passes):
        for k in range(K):
            def res(p):
                bb = b.copy(); bb[k] = np.exp(p[:4])
                aa = a.copy(); aa[k] = np.exp(p[4])
                z = sim(bb, rho, a=aa)
                return np.sqrt(z[:, k]) - np.sqrt(ys[:, k])
            starts = [np.r_[np.log(b[k]), np.log(a[k])]] + [np.r_[np.log(np.full(4, s)), 0.0] for s in (0.004, 0.01, 0.025)]
            best = None
            lo = np.r_[np.log(np.full(4, 5e-4)), np.log(0.05)]
            hi = np.r_[np.log(np.full(4, 0.5)), np.log(20)]
            for p0 in starts:
                r = least_squares(res, p0, bounds=(lo, hi))
                if best is None or r.cost < best.cost:
                    best = r
            b[k] = np.exp(best.x[:4]); a[k] = np.exp(best.x[4])
        if verbose:
            print("passagem", it, "custo", round(cost(sim(b, rho, a=a)), 2), flush=True)
    return b, a


def metrics(z):
    rmse = np.sqrt(((z - ys) ** 2).mean(0))
    r2 = 1 - ((z - ys) ** 2).sum(0) / ((ys - ys.mean(0)) ** 2).sum(0)
    return rmse, r2


def R0_block(k, bk):
    """R0 do bloco regional (eq. 5.6): T_I * raio espectral de F_kk."""
    M = np.array([[N[k, g] * bk * Th[k, k] ** 2 * C[g, h] / Nhat[k] for h in range(G)] for g in range(G)])
    return TI * max(abs(np.linalg.eigvals(M)))


def load_params(rho):
    return np.load(CALIB / f"beta_rho{rho}.npy"), np.load(CALIB / f"a_rho{rho}.npy")
