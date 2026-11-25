import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from seirmoz.control import setup, K, G


def test_obitos_sem_controlo_e_com_controlo():
    P, J0, k, b = setup(0.1)
    z = np.zeros((P.n, K, G))
    X0, _ = P.forward(z, z)
    d0 = (X0[0][10] + X0[0][11]).sum(); dT = (X0[-1][10] + X0[-1][11]).sum()
    assert abs((dT - d0) - 17673) < 5
    v, u, X, J, h, it = P.solve()
    d = (X[-1][10] + X[-1][11]).sum() - d0
    assert abs(d - 3921) < 5
