"""Figures for the profiles document."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.colors import PowerNorm
import toy
from core import axes, particle_profiles, mode_profiles

plt.rcParams.update({"font.size": 8.5, "mathtext.fontset": "cm", "font.family": "serif",
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 9, "axes.titlelocation": "left"})
BLUE, AQUA, ORANGE, GRAY, VIOLET, INK = "#2a78d6", "#1baf7a", "#eb6834", "#8a8984", "#6250d6", "#52514e"
r, N, M, springs = toy.r, toy.N, toy.M, toy.springs
i_idx, n_idx = np.arange(1, N + 1), np.arange(1, M + 1)


# ------------------------------------------------------------------ worked examples (section 4)
def draw_ellipsoid(ax, C, color):
    s2, u = axes(C)
    s = np.sqrt(s2)
    a, b = np.meshgrid(np.linspace(0, 2 * np.pi, 40), np.linspace(0, np.pi, 20))
    pts = (np.stack([np.cos(a) * np.sin(b), np.sin(a) * np.sin(b), np.cos(b)], axis=-1) * s) @ u.T
    ax.plot_wireframe(pts[..., 0], pts[..., 1], pts[..., 2], color=color, linewidth=0.35, alpha=0.6)
    for k in range(3):
        v = u[:, k] * s[k]
        ax.plot(*np.stack([-v, v]).T, color=color, linewidth=1.6)


def rot_z(deg):
    t = np.radians(deg)
    return np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])


def examples():
    needle = np.diag(np.array([4 ** (2 / 3), 4 ** (-1 / 3), 4 ** (-1 / 3)]) ** 2)
    cases = [("(a) size", np.eye(3), 4 * np.eye(3)),
             ("(b) shape", np.eye(3), needle),
             ("(c) orientation", needle, rot_z(90) @ needle @ rot_z(90).T)]
    fig = plt.figure(figsize=(7.2, 4.9))
    for k, (title, CA, CB) in enumerate(cases):
        ax = fig.add_axes([0.0 + 0.33 * k, 0.46, 0.33, 0.44], projection="3d")
        draw_ellipsoid(ax, CA, BLUE)
        draw_ellipsoid(ax, CB, AQUA)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlim(-2.6, 2.6); ax.set_ylim(-2.6, 2.6); ax.set_zlim(-2.6, 2.6)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.view_init(elev=22, azim=-55)
        p = particle_profiles(CA, CB)
        ax.set_title(f"{title}\n$d_{{\\mathrm{{size}}}}={p['d_size'][0]:.2f}$,  $d_{{\\mathrm{{shape}}}}={p['d_shape'][0]:.2f}$,"
                     f"  $d_{{\\mathrm{{orient}}}}={p['d_orient'][0]:.2f}$", fontsize=8.5, loc="center")
    fig.legend(handles=[plt.Line2D([], [], color=BLUE, lw=2, label="A"), plt.Line2D([], [], color=AQUA, lw=2, label="B")],
               loc="upper center", bbox_to_anchor=(0.5, 0.49), ncol=2, frameon=False)
    ax = fig.add_axes([0.12, 0.09, 0.5, 0.3])
    th = np.linspace(0, 90, 181)
    for q, shade in zip([1.1, 1.5, 2, 4], ["#a9c9f0", "#6fa3e3", BLUE, "#13427a"]):
        CA = np.diag([q ** 2, 1, 1])
        d = [particle_profiles(CA, rot_z(t) @ CA @ rot_z(t).T)["d_orient"][0] for t in th]
        ax.plot(th, d, color=shade, lw=1.6, label=f"{q:g}:1:1")
    ax.set_xlabel("rotation angle of B relative to A (degrees)")
    ax.set_ylabel(r"$d_{\mathrm{orient}}$")
    ax.set_xlim(0, 90); ax.set_ylim(0, None); ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
    ax.set_title("(d) rotation of a needle about a short axis", fontsize=8.5)
    ax.legend(title="standard deviations", frameon=False, fontsize=7.5, title_fontsize=7.5,
              loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.savefig("fig_examples.pdf")
    plt.close(fig)


# ------------------------------------------------------------------ toy model and reference system (section 6)
def network(ax, mask, lw_other=0.4):
    for (i, j), hit in zip(springs, mask):
        if not hit:
            ax.plot(r[[i, j], 0], r[[i, j], 1], color="#dcdbd6", lw=lw_other, zorder=1)
    for (i, j), hit in zip(springs, mask):
        if hit:
            ax.plot(r[[i, j], 0], r[[i, j], 1], color=ORANGE, lw=0.8, zorder=2)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def domain_line(ax):
    ax.axvline(N // 2 + 0.5, color=GRAY, lw=0.6, ls="--")
    ax.set_xlim(0.5, N + 0.5)


def toy_figure():
    s2A = toy.sigma2_A
    p = particle_profiles(toy.C_A, toy.C_A)
    fig = plt.figure(figsize=(7.4, 6.0))
    for k, (key, cfg) in enumerate(toy.PERTURBATIONS.items()):
        ax = fig.add_axes([0.02 + 0.33 * k, 0.70, 0.3, 0.26])
        network(ax, cfg["mask"])
        ax.scatter(r[:, 0], r[:, 1], s=8, c=np.where(toy.domain1, "#3b3a37", "#9b9a95"), zorder=3, linewidths=0)
        ax.set_title(f"({'abc'[k]}) {cfg['name']}", fontsize=8.5)
    rows = [("mu_A", r"$\mu^A_i$", "(d) size"), ("e_A", r"$e^A_i$", "(e) elongation"), ("f_A", r"$f^A_i$", "(f) flatness")]
    for k, (col, ylab, title) in enumerate(rows):
        ax = fig.add_axes([0.08, 0.53 - 0.205 * k, 0.5, 0.14])
        ax.plot(i_idx, p[col], "-o", color=BLUE, ms=2.5, lw=1.2)
        domain_line(ax)
        ax.set_ylabel(ylab); ax.set_title(title, fontsize=8.5)
        ax.grid(alpha=0.25, lw=0.5)
        if k == 2:
            ax.set_xlabel(r"particle $i$")
    ax = fig.add_axes([0.7, 0.11, 0.27, 0.56])
    ax.plot(n_idx, np.log(s2A), color=BLUE, lw=1.6)
    ax.set_xlabel(r"mode $n$"); ax.set_ylabel(r"$\ln(\sigma^A_n)^2$")
    ax.set_title("(g) spectrum", fontsize=8.5)
    ax.set_xlim(0, M + 1)
    ax.grid(alpha=0.25, lw=0.5)
    fig.savefig("fig_toy.pdf")
    plt.close(fig)


# ------------------------------------------------------------------ one page per perturbation (section 7)
def signed_bars(ax, x, v):
    ax.bar(x, v, width=0.85, color=np.where(v >= 0, AQUA, BLUE), linewidth=0)
    ax.axhline(0, color=GRAY, lw=0.6)
    pad = 0.1 * np.max(np.abs(v))
    ax.set_ylim(min(0.0, v.min()) - pad, max(0.0, v.max()) + pad)


def mark(ax, key):
    if key == "interface":
        cfg = toy.PERTURBATIONS[key]
        ii = np.unique(springs[cfg["mask"]]) + 1
        ax.plot(ii, np.zeros(len(ii)), "|", color=ORANGE, ms=6, mew=1.2,
                transform=ax.get_xaxis_transform(), clip_on=False)
    else:
        ax.axvspan(8.5, 13.5, color="#f8e1d6", zorder=0)


def ellipses(ax, C, color, scale):
    for i in range(N):
        b = slice(3 * i, 3 * i + 3)
        w, V = np.linalg.eigh(C[b, b][:2, :2])
        ang = np.degrees(np.arctan2(V[1, 1], V[0, 1]))
        ax.add_patch(Ellipse(r[i, :2], 2 * scale * np.sqrt(w[1]), 2 * scale * np.sqrt(w[0]),
                             angle=ang, fill=False, color=color, lw=0.9, zorder=4))


def perturbation_figure(key):
    cfg = toy.PERTURBATIONS[key]
    s2B, UB, CB = toy.system_B(key)
    p = particle_profiles(toy.C_A, CB)
    m = mode_profiles(toy.sigma2_A, toy.U_A, CB, UB)
    fig = plt.figure(figsize=(7.4, 9.6))

    ax = fig.add_axes([0.1, 0.745, 0.87, 0.235])
    network(ax, cfg["mask"], lw_other=0.3)
    ellipses(ax, toy.C_A, BLUE, 2.5)
    ellipses(ax, CB, AQUA, 2.5)
    ax.autoscale_view()
    ax.set_title(r"(a) particle covariances $C_{ii}$", fontsize=8.5)

    rows = [(p["mu_B"] - p["mu_A"], r"$\Delta\mu_i$", "(b) size"),
            (p["e_B"] - p["e_A"], r"$\Delta e_i$", "(c) elongation"),
            (p["f_B"] - p["f_A"], r"$\Delta f_i$", "(d) flatness"),
            (p["d_orient"], r"$d_{\mathrm{orient},i}$", "(e) orientation")]
    for k, (v, ylab, title) in enumerate(rows):
        ax = fig.add_axes([0.1, 0.64 - 0.093 * k, 0.87, 0.064])
        mark(ax, key)
        if k < 3:
            signed_bars(ax, i_idx, v)
        else:
            ax.bar(i_idx, v, width=0.85, color=GRAY, linewidth=0)
            ax.set_ylim(0, 1.1 * v.max())
        domain_line(ax)
        ax.set_ylabel(ylab); ax.set_title(title, fontsize=8.5)
        if k == 3:
            ax.set_xlabel(r"particle $i$")
        else:
            ax.set_xticklabels([])

    ax = fig.add_axes([0.1, 0.205, 0.5, 0.1])
    ax.bar(n_idx, m["N"], width=0.85, color=VIOLET, linewidth=0)
    ax.axhline(1, color=GRAY, lw=0.6)
    ax.set_xlim(0.5, M + 0.5); ax.set_xticklabels([])
    ax.set_ylabel(r"$N_n$"); ax.set_title("(f) conservation of each mode of A", fontsize=8.5)
    ax = fig.add_axes([0.1, 0.055, 0.5, 0.1])
    signed_bars(ax, n_idx, m["rho"])
    ax.set_xlim(0.5, M + 0.5)
    ax.set_xlabel(r"mode $n$ of A"); ax.set_ylabel(r"$\rho_n$")
    ax.set_title("(g) fluctuation change along each mode of A", fontsize=8.5)

    ax = fig.add_axes([0.70, 0.055, 0.235, 0.25])
    im = ax.imshow(m["O2"], cmap="Purples", norm=PowerNorm(0.5, vmin=0, vmax=1), interpolation="nearest",
                   extent=[0.5, M + 0.5, M + 0.5, 0.5])
    ax.set_xlabel(r"mode $m$ of B"); ax.set_ylabel(r"mode $n$ of A")
    ax.set_title(r"(h) $O^2_{nm}$", fontsize=8.5)
    ax.set_xticks([1, 50, 100, 144]); ax.set_yticks([1, 50, 100, 144])
    cax = fig.add_axes([0.95, 0.08, 0.008, 0.2])
    cb = fig.colorbar(im, cax=cax, ticks=[0, 0.05, 0.2, 0.5, 1]); cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7)
    fig.savefig(f"fig_{key}.pdf")
    fig.savefig(f"fig_{key}.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    examples()
    toy_figure()
    for key in toy.PERTURBATIONS:
        perturbation_figure(key)
