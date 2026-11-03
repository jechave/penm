"""Per-particle and per-mode comparison of two elastic network models."""
import numpy as np


def hessian(r, springs, k):
    """3N x 3N Hessian of an elastic network: positions r (N x 3), springs (pairs i, j), constants k."""
    N = len(r)
    H = np.zeros((3 * N, 3 * N))
    for (i, j), kij in zip(springs, k):
        d = r[j] - r[i]
        b = kij * np.outer(d, d) / (d @ d)
        H[3*i:3*i+3, 3*j:3*j+3] -= b
        H[3*j:3*j+3, 3*i:3*i+3] -= b
        H[3*i:3*i+3, 3*i:3*i+3] += b
        H[3*j:3*j+3, 3*j:3*j+3] += b
    return H


def modes(H, n_rigid=6):
    """Variances sigma_n^2 (softest first), modes U (3N x M) and covariance C = pseudo-inverse of H."""
    w, V = np.linalg.eigh(H)                      # ascending: rigid-body modes first, then softest to stiffest
    sigma2, U = 1.0 / w[n_rigid:], V[:, n_rigid:]
    return sigma2, U, (U * sigma2) @ U.T


def axes(C_ii):
    """Variances sigma_{i,alpha}^2 (longest axis first) and axes u_{i,alpha} of a 3 x 3 covariance."""
    s2, u = np.linalg.eigh(C_ii)
    return s2[::-1], u[:, ::-1]


def particle_profiles(C_A, C_B):
    """For every particle i: mu, e, f in A and in B, and d_size, d_shape, d_orient."""
    out = []
    for i in range(C_A.shape[0] // 3):
        blk = slice(3 * i, 3 * i + 3)
        s2_A, u_A = axes(C_A[blk, blk])
        s2_B, u_B = axes(C_B[blk, blk])
        mu_A, mu_B = np.log(s2_A).mean(), np.log(s2_B).mean()
        l_A, l_B = np.log(s2_A) - mu_A, np.log(s2_B) - mu_B
        o2 = (u_A.T @ u_B) ** 2                    # squared overlaps o_{i,alpha beta}^2
        d_size2 = 3 * (mu_B - mu_A) ** 2
        d_shape2 = np.sum((l_B - l_A) ** 2)
        d_orient2 = max(0.0, 2 * l_A @ (np.eye(3) - o2) @ l_B)
        out.append([mu_A, l_A[0] - l_A[1], l_A[1] - l_A[2],
                    mu_B, l_B[0] - l_B[1], l_B[1] - l_B[2],
                    np.sqrt(d_size2), np.sqrt(d_shape2), np.sqrt(d_orient2)])
    names = ["mu_A", "e_A", "f_A", "mu_B", "e_B", "f_B", "d_size", "d_shape", "d_orient"]
    return dict(zip(names, np.array(out).T))


def mode_profiles(sigma2_A, U_A, C_B, U_B):
    """For every mode n of A: N_n and rho_n; also the squared overlaps O_nm^2."""
    O2 = (U_A.T @ U_B) ** 2
    S = -np.sum(np.where(O2 > 0, O2 * np.log(np.where(O2 > 0, O2, 1.0)), 0.0), axis=1)
    rho = np.log(np.einsum("in,ij,jn->n", U_A, C_B, U_A)) - np.log(sigma2_A)
    return dict(N=np.exp(S), rho=rho, O2=O2)
