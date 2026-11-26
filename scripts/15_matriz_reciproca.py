"""Extensão: matriz de contactos simetrizada (N_g C_gh = N_h C_hg), recalibração e controlo.
Uso: python scripts/15_matriz_reciproca.py cal | ctl"""
import sys, pathlib, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from seirmoz import calibration as cal, model as M
from seirmoz.paths import CALIB, RESULTS
Ng = M.N.sum(0); Cs = (Ng[:, None] * M.C + (Ng[:, None] * M.C).T) / (2 * Ng[:, None])
cal.C[:] = Cs; M.C[:] = Cs          # as duas cópias são o mesmo objecto em ambos os módulos
D = CALIB.parent / 'calibracao_recip'; D.mkdir(exist_ok=True)
cal.calib_dir = lambda: D
if sys.argv[1] == 'cal':
    for rho in (0.1,):
        t0 = time.time(); beta, a = cal.fit(rho, passes=3)
        np.save(D / f"beta_rho{rho}.npy", beta); np.save(D / f"a_rho{rho}.npy", a)
        rmse, r2 = cal.metrics(cal.sim(beta, rho, a))
        print("rho", rho, "tempo", round(time.time() - t0), "s\nbeta\n", np.round(beta, 4), "\nR2", np.round(r2, 3), "\nRMSE", np.round(rmse, 1), flush=True)
else:
    from seirmoz.control import setup, K, G, N
    out = {}
    for nome, mode in (("recip", 1),):
        P, J0, k, b = setup(0.1); z = np.zeros((P.n, K, G)); X0, _ = P.forward(z, z)
        v, u, X, J, h, it = P.solve(it=400); d0 = P.x0[10] + P.x0[11]
        ob0 = (X0[-1][10] + X0[-1][11] - d0).sum(); ob = (X[-1][10] + X[-1][11] - d0).sum()
        print(nome, "kappa", round(k, 3), "sem", round(ob0), "com", round(ob), "red", round(100 * (ob0 - ob) / ob0, 1), "doses", round((X[-1][12] - P.x0[12]).sum()),
              "cob idade", np.round(100 * (X[-1][12] - P.x0[12]).sum(0) / N.sum(0), 1), "u idade", np.round((u.mean(0) * N).sum(0) / N.sum(0), 3), flush=True)
