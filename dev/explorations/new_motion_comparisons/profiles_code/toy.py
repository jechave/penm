"""Toy two-domain elastic network and its perturbations."""
import numpy as np
from core import hessian, modes


def grow_chain(n, center, radius, rng, start, avoid=None, bond=3.8, min_dist=4.0):
    avoid = np.zeros((0, 3)) if avoid is None else avoid
    for _ in range(2000):
        r = [np.array(start, float)]
        for _ in range(n - 1):
            for _ in range(500):
                d = rng.normal(size=3)
                p = r[-1] + bond * d / np.linalg.norm(d)
                if np.linalg.norm(p - center) > radius:
                    continue
                others = np.vstack([avoid] + r[:-1]) if (len(avoid) or len(r) > 1) else np.zeros((0, 3))
                if len(others) and np.min(np.linalg.norm(others - p, axis=1)) < min_dist:
                    continue
                r.append(p)
                break
            else:
                break
        if len(r) == n:
            return np.array(r)
    raise RuntimeError("chain growth failed")


def two_domain_structure(n_per_domain=25, separation=15.0, radius=7.5, seed=3):
    rng = np.random.default_rng(seed)
    c1, c2 = np.array([-separation / 2, 0, 0]), np.array([separation / 2, 0, 0])
    r1 = grow_chain(n_per_domain, c1, radius, rng, start=c1)
    r2 = grow_chain(n_per_domain, c2, radius, rng, start=c2, avoid=r1)
    return np.vstack([r1, r2])


def springs_within(r, R_c):
    N = len(r)
    return np.array([(i, j) for i in range(N) for j in range(i + 1, N) if np.linalg.norm(r[i] - r[j]) < R_c])


r = two_domain_structure()
N = len(r)
R_C = 10.0
springs = springs_within(r, R_C)
k_A = np.ones(len(springs))
sigma2_A, U_A, C_A = modes(hessian(r, springs, k_A))
M = len(sigma2_A)

domain1 = np.arange(N) < N // 2
interface = domain1[springs[:, 0]] != domain1[springs[:, 1]]
site = np.arange(8, 13)                                   # particles 9-13
local = np.isin(springs[:, 0], site) | np.isin(springs[:, 1], site)
z = np.random.default_rng(11).normal(size=len(springs))

GAMMA, ETA = 0.3, 0.8
PERTURBATIONS = {
    "local": dict(name="local softening", mask=local, factor=np.where(local, GAMMA, 1.0)),
    "interface": dict(name="interface softening", mask=interface, factor=np.where(interface, GAMMA, 1.0)),
    "random": dict(name="random local springs", mask=local, factor=np.where(local, np.exp(ETA * z), 1.0)),
}


def system_B(key):
    return modes(hessian(r, springs, k_A * PERTURBATIONS[key]["factor"]))
