"""Extensão: calibração com a convenção padrão C_gh I_h / N_h (denominador N̂_{mh})."""
import sys, pathlib, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from seirmoz import calibration as cal, model

cal.set_norm('idade'); model.set_norm('idade')
rhos = [float(x) for x in sys.argv[1:]] or [0.1]
for rho in rhos:
    t0 = time.time()
    beta, a = cal.fit(rho, passes=3)
    d = cal.calib_dir(); d.mkdir(exist_ok=True)
    np.save(d / f"beta_rho{rho}.npy", beta); np.save(d / f"a_rho{rho}.npy", a)
    z = cal.sim(beta, rho, a)
    rmse, r2 = cal.metrics(z)
    np.set_printoptions(precision=4, suppress=True)
    print("rho", rho, "tempo", round(time.time() - t0), "s")
    print("beta\n", beta, "\nR2", r2, "\nRMSE", rmse, flush=True)
