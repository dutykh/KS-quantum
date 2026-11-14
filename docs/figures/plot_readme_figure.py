#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: quasinormal modes of the Kazakov-Solodukhin black hole read from reports/.

From this folder,
    python3 plot_readme_figure.py
writes ks-quantum-modes.png. Panels (a) and (b): the fundamental mode of four sectors against the
horizon-normalised deformation a_K = a / r_h. Panel (c): every retained root of the scalar monopole
at a_K = 0.9. The style is the stamped copy manuscript_style.py beside this script.
"""
import math
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPORTS = HERE.parent.parent / "reports"
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

SECTORS = (  # folder, file pattern, legend
    ("s0", "s0_L0_*.txt", r"scalar, $\ell = 0$"),
    ("s0-conformal", "s0_L0_*.txt", r"conformal scalar, $\ell = 0$"),
    ("s1", "s1_L1_*.txt", r"electromagnetic, $\ell = 1$"),
    ("s2-axial", "s2_L2_*.txt", r"axial, $\ell = 2$"),
)
SMALLEST = 1e-3   # roots of smaller modulus are left out of panel (c)


def read_report(path):
    """Deformation a = a/M from the header, and the retained roots of the table."""
    text = path.read_text()
    a = float(re.search(r"a-tag = ([0-9.]+)", text).group(1)) / 1000
    roots = []
    for line in text.splitlines():
        p = line.split()
        if len(p) == 7 and p[0].isdigit():
            roots.append(complex(float(p[1]), float(p[2])))
    return a, np.array(roots)


def fundamental_curve(folder, pattern):
    """(a_K, Omega) of the least damped root with Re Omega > 0, sorted by a_K."""
    points = {}
    for path in sorted((REPORTS / folder).glob(pattern)):
        a, roots = read_report(path)
        oscillatory = roots[roots.real > SMALLEST]
        if oscillatory.size == 0:
            continue
        a_K = round(a / math.sqrt(a * a + 4.0), 4)
        points.setdefault(a_K, oscillatory[np.argmax(oscillatory.imag)])
    a_K = np.array(sorted(points))
    omega = np.array([points[k] for k in a_K])
    # A gap wider than 0.2 in a_K is left open rather than bridged by a straight line.
    gaps = np.flatnonzero(np.diff(a_K) > 0.2)
    return np.insert(a_K, gaps + 1, np.nan), np.insert(omega, gaps + 1, np.nan)


hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 265
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig, (ax_re, ax_im, ax_sp) = plt.subplots(1, 3, figsize=(hs.WIDTH, 2.45),
                                          gridspec_kw={"width_ratios": [1, 1, 1.05]})
for k, (folder, pattern, label) in enumerate(SECTORS):
    a_K, omega = fundamental_curve(folder, pattern)
    style = hs.series(k)
    style["markersize"] = 2.6
    ax_re.plot(a_K, omega.real, lw=0.9, label=label, **style)
    ax_im.plot(a_K, -omega.imag, lw=0.9, **style)

for axis in (ax_re, ax_im):
    axis.set_xlim(0.0, 1.0)
    axis.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0], ["0", "0.25", "0.5", "0.75", "1"])
    axis.set_xlabel(r"$a_K = a / r_h$")
ax_re.set_ylim(0.0, 0.4)
ax_re.set_ylabel(r"$\mathrm{Re}\,\Omega$")
hs.panel(ax_re, "(a) frequency")
ax_im.set_ylim(0.0, 0.12)
ax_im.set_ylabel(r"$-\mathrm{Im}\,\Omega$")
hs.panel(ax_im, "(b) damping rate")

_, roots = read_report(REPORTS / "s0" / "s0_L0_aK0_9.txt")
roots = roots[np.abs(roots) > SMALLEST]
ax_sp.plot(roots.real, roots.imag, linestyle="none", marker="o", markersize=2.6,
           markerfacecolor="white", markeredgecolor=hs.NAVY, markeredgewidth=0.8)
ax_sp.set_xlim(-0.14, 0.14)
ax_sp.set_ylim(-15.5, 0.0)
ax_sp.set_xticks([-0.1, 0.0, 0.1])
ax_sp.set_xlabel(r"$\mathrm{Re}\,\Omega$")
ax_sp.set_ylabel(r"$\mathrm{Im}\,\Omega$")
hs.panel(ax_sp, r"(c) scalar, $\ell = 0$, $a_K = 0.9$")

hs.legend_outside(fig, where="upper", ncol=4)
out = HERE / "ks-quantum-modes.png"
hs.save(fig, out)
print("written", out)
