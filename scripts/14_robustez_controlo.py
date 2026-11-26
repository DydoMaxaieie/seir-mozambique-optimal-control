"""Extensão: (a) condições iniciais do FBS; (b) sensibilidade à eficácia vacinal e aos parâmetros clínicos."""
import sys, pathlib, csv, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from seirmoz import model as M
from seirmoz.control import setup, K, G, N
from seirmoz.paths import RESULTS

def clinico(ihr_scale=1.0, phi_out=0.30, eff=0.70, ifr_exp=1.0):
    """Recalcula os parâmetros clínicos derivados (mesmas fórmulas de model.py)."""
    M.EFF = eff; M.PHI_OUT = phi_out
    ihr = M.IHR * ihr_scale
    ph = np.maximum(M.Wt @ ihr, 1e-4); ifr = M.Wt @ (M.IFR ** ifr_exp * 0.006 ** (1 - ifr_exp))
    pc = (M.Wt @ (ihr * M.PC)) / np.maximum(M.Wt @ ihr, 1e-12); pc[0] = 0.05
    fq = phi_out * ifr / (1 - ph); h_ = (1 - phi_out) * ifr / ph
    fw = 0.3 * h_; fc = np.minimum((h_ - fw) / ((1 - fw) * pc), 0.5)
    M.ph, M.ifr, M.pc, M.fq, M.h_, M.fw, M.fc = ph, ifr, pc, fq, h_, fw, fc
    M.mu_q = fq / ((1 - fq) * M.TQ); M.mu_w = fw / ((1 - fw) * M.THW); M.mu_c = fc / ((1 - fc) * M.THC)

def corre(rho=0.1, label=""):
    P, J0, k, b = setup(rho)
    z = np.zeros((P.n, K, G)); X0, Jz = P.forward(z, z)
    v, u, X, J, h, it = P.solve(it=400)
    d0 = P.x0[10] + P.x0[11]
    ob0 = (X0[-1][10] + X0[-1][11] - d0).sum(); ob = (X[-1][10] + X[-1][11] - d0).sum()
    return dict(label=label, kappa=k, sem=ob0, vu=ob, red=100 * (ob0 - ob) / ob0, doses=(X[-1][12] - P.x0[12]).sum(), J_J0=J / J0, it=it,
                cob=((X[-1][12] - P.x0[12]).sum(0) / N.sum(0)), v=v, u=u, P=P)

rows = []
modo = sys.argv[1] if len(sys.argv) > 1 else "ic"
if modo == "ic":
    base = corre(label="base (0)")
    P = base["P"]; rng = np.random.default_rng(1)
    print(f"base: J/J0 {base['J_J0']:.5f} óbitos {base['vu']:.1f} it {base['it']}", flush=True)
    inits = {"v=vmax/2, u=umax/2": (np.broadcast_to(P.vmax / 2, (P.n, K, G)).copy(), np.full((P.n, K, G), P.umax / 2)),
             "v=vmax, u=umax": (np.broadcast_to(P.vmax, (P.n, K, G)).copy(), np.full((P.n, K, G), P.umax)),
             "aleatório 1": (rng.uniform(0, 1, (P.n, K, G)) * P.vmax, rng.uniform(0, P.umax, (P.n, K, G))),
             "aleatório 2": (rng.uniform(0, 1, (P.n, K, G)) * P.vmax, rng.uniform(0, P.umax, (P.n, K, G)))}
    for nome, (v0, u0) in inits.items():
        v, u, X, J, h, it = P.solve(v0=v0, u0=u0, it=400)
        d0 = P.x0[10] + P.x0[11]; ob = (X[-1][10] + X[-1][11] - d0).sum()
        dv = np.abs(v - base["v"]).max() / P.vmax.max(); du = np.abs(u - base["u"]).max() / P.umax
        print(f"{nome:22s} J/J0 {J/ (base['J_J0']**-1 * 0 + 1) :.5f} (rel. base {J/ (base['J_J0']*1):.0f}) óbitos {ob:.1f} it {it} max|dv|/vmax {dv:.4f} max|du|/umax {du:.4f}", flush=True)
        rows.append([nome, J, ob, it, dv, du])
    with open(RESULTS / "controlo" / "condicoes_iniciais.csv", "w", newline="", encoding="utf8") as f:
        w = csv.writer(f); w.writerow(["inicial", "J", "obitos", "iteracoes", "dv_rel", "du_rel"]); w.writerows(rows)
else:
    casos = [("base", {}), ("eficácia vacinal 0,5", dict(eff=0.5)), ("eficácia vacinal 0,9", dict(eff=0.9)),
             ("IHR x0,5", dict(ihr_scale=0.5)), ("IHR x2", dict(ihr_scale=2.0)),
             ("gradiente etário da IFR mais plano (expoente 0,8)", dict(ifr_exp=0.8)), ("gradiente etário da IFR mais íngreme (expoente 1,2)", dict(ifr_exp=1.2)),
             ("óbitos fora do hospital 10%", dict(phi_out=0.10)), ("óbitos fora do hospital 50%", dict(phi_out=0.50))]
    for nome, kw in casos:
        clinico(**kw); t = time.time(); r = corre(label=nome)
        print(f"{nome:30s} kappa {r['kappa']:.2f} sem {r['sem']:8.0f} com {r['vu']:7.0f} red {r['red']:5.1f}% doses {r['doses']/1e3:6.0f}k cob/idade {np.round(100*r['cob'][2:],1)} it {r['it']} ({time.time()-t:.0f}s)", flush=True)
        rows.append([nome, round(r['kappa'], 3), round(r['sem'], 1), round(r['vu'], 1), round(r['red'], 1), round(r['doses']), *np.round(100 * r['cob'][2:], 1)])
    with open(RESULTS / "controlo" / "sensibilidade_clinica.csv", "w", newline="", encoding="utf8") as f:
        w = csv.writer(f); w.writerow(["caso", "kappa", "obitos_sem", "obitos_vu", "reducao_pct", "doses", "cob_18_39", "cob_40_59", "cob_60"]); w.writerows(rows)
