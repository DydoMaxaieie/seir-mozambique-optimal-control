"""Entradas demográficas, de contacto e de mobilidade (Secção 5.1 da dissertação).

Constrói N_kg (população), C (contactos agregados, eq. 5.1) e Theta (mobilidade, eq. 5.2)
a partir dos ficheiros de data/.
"""
import csv
from functools import lru_cache
import numpy as np
from .paths import DATA

REGS = ["Norte", "Centro", "Sul-interior", "Sul-metropolitano"]
MR = {
    "Norte": ["Niassa", "Cabo Delgado", "Nampula", "Zambezia"],
    "Centro": ["Tete", "Manica", "Sofala"],
    "Sul-interior": ["Inhambane", "Gaza"],
    "Sul-metropolitano": ["Maputo", "Maputo Cidade"],
}
# nomes das províncias na matriz de fluxos (Joaquim et al., 2025) -> chave de MR
_FLUX_ORDER = ["Niassa", "Cabo Delgado", "Nampula", "Zambezia", "Tete", "Manica", "Sofala",
               "Inhambane", "Gaza", "Maputo", "Maputo Cidade"]


def _read_prov17():
    with open(DATA / "populacao_provincias_2017.csv", encoding="utf8") as f:
        return {r["provincia"]: float(r["populacao_2017"]) for r in csv.DictReader(f)}


def _read_flows():
    with open(DATA / "fluxos_provinciais.csv", encoding="utf8") as f:
        rows = list(csv.reader(f))[1:]
    return np.array([[float(x) for x in r[1:]] for r in rows])


@lru_cache(maxsize=1)
def load():
    prov17 = _read_prov17()
    # população por classe de 5 anos (ONU 2019, WPP), 0-4 ... 75+ (16 classes)
    with open(DATA / "populacao_idade_2020.csv") as f:
        un = np.array([float(r["populacao_milhares"]) for r in csv.DictReader(f)]) * 1000
    bins = np.r_[un[:15], un[15:].sum()]
    tot = bins.sum()
    share17 = sum(prov17.values())
    Nk = np.array([sum(prov17[p] for p in MR[r]) / share17 * tot for r in REGS])
    # grupos 0-4, 5-17, 18-39, 40-59, 60+ (classe 15-19 repartida 3/5 e 2/5)
    frac = np.zeros((16, 5))
    frac[0, 0] = 1; frac[1, 1] = 1; frac[2, 1] = 1; frac[3, 1] = .6; frac[3, 2] = .4
    for i in range(4, 8): frac[i, 2] = 1
    for i in range(8, 12): frac[i, 3] = 1
    for i in range(12, 16): frac[i, 4] = 1
    Ng = bins @ frac
    ag = Ng / Ng.sum()
    N = np.outer(Nk, ag)
    # contactos (Prem et al., 2021), agregados por (5.1)
    c = np.loadtxt(DATA / "matriz_contactos_16x16.csv", delimiter=",")
    W = bins[:, None] * frac
    C = np.zeros((5, 5))
    for g in range(5):
        for h in range(5):
            C[g, h] = sum(W[i, g] * (c[i, :] @ frac[:, h]) for i in range(16)) / W[:, g].sum()
    # mobilidade: soma dos fluxos entre macrorregiões e eq. (5.2)
    T = _read_flows()
    idx = {r: [_FLUX_ORDER.index(p) for p in MR[r]] for r in REGS}
    A = np.zeros((4, 4))
    for a, ra in enumerate(REGS):
        for b, rb in enumerate(REGS):
            if a != b:
                A[a, b] = T[np.ix_(idx[ra], idx[rb])].sum()
    Th = A / Nk[:, None]
    for k in range(4):
        Th[k, k] = 1 - Th[k].sum()
    return dict(N=N, C=C, Th=Th, Nk=Nk, ag=ag, A=A)


if __name__ == "__main__":
    d = load()
    np.set_printoptions(suppress=True, linewidth=120)
    print("Nk (milhões):", np.round(d["Nk"] / 1e6, 2), "total", round(d["Nk"].sum() / 1e6, 2))
    print("proporção por grupo:", np.round(d["ag"], 3))
    print("C =\n", np.round(d["C"], 2))
    print("Theta =\n", np.round(d["Th"], 5))
