# Quasinormal modes of the Kazakov–Solodukhin quantum-corrected black hole — a spectral analysis

This repository contains the data, symbolic worksheets, and numerical code that accompany
the manuscript

> **Quasinormal modes of the Kazakov–Solodukhin quantum-corrected black hole: a spectral analysis**
> Davide Batic, Denys Dutykh, and Mark Essa Sukaiti
> *Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, United Arab Emirates*

A copy of the paper is included here: [`pdf/QNM_Kazakov_Solodukhin_quantum_corrected_black_hole.pdf`](pdf/QNM_Kazakov_Solodukhin_quantum_corrected_black_hole.pdf).

We compute the quasinormal-mode (QNM) spectrum of the **Kazakov–Solodukhin quantum-corrected
(KSQC) black hole** with a high-precision **Chebyshev spectral method**. After factoring out the
QNM asymptotics at the event horizon and at spatial infinity, the radial problem reduces to a
**quadratic matrix pencil** for the dimensionless frequency $\Omega = M\omega$. We treat four
perturbation sectors — minimally coupled scalar, electromagnetic, non-minimally (conformally)
coupled scalar, and an axial effective-source Regge–Wheeler gravitational model — and validate
the method against the WKB, time-domain, Mashhoon, and asymptotic-iteration results in the
literature.

> **Sign convention.** Throughout, the time dependence is $e^{-i\omega t}$, so a damped
> (physical) mode has $\mathrm{Im}\,\Omega < 0$. A positive imaginary part would signal an
> instability; within the parameter ranges and resolutions explored here, **no growing mode is
> found**.

---

## 1. The physical system

The KSQC black hole is a spherically symmetric, **non-Ricci-flat** deformation of the
Schwarzschild solution. In Planck units ($c = G_N = \hbar = 1$) the line element is

$$ds^2 = -f(r)\,dt^2 + f(r)^{-1}\,dr^2 + r^2\,(d\vartheta^2 + \sin^2\vartheta\,d\varphi^2),
\qquad f(r) = \frac{1}{r}\Big(\sqrt{r^2 - a^2} - 2M\Big), \quad r > a,$$

where $M$ is the mass and $a^2 = 4 G_R/\pi$ is the **deformation parameter** measuring the
strength of the underlying quantum correction. The event horizon sits at
$r_h = \sqrt{4M^2 + a^2}$, and the geometry carries a curvature singularity at $r = a$. Because
the spacetime is not Ricci flat, the gravitational sector requires special care (see below).

The radial perturbation equations for the various sectors are reduced to a single dimensionless
form for $\Omega = M\omega$ and solved spectrally. The four sectors treated are:

| Sector | Spin $s$ | Field | Notes |
|---|---|---|---|
| Minimally coupled scalar | 0 | massless scalar | $\ell \ge 0$ |
| Electromagnetic | 1 | vector | $\ell \ge 1$ |
| Non-minimally coupled scalar | 0 | conformally coupled scalar | adds the $R/6$ curvature term to the potential |
| Axial gravitational | 2 | effective-source Regge–Wheeler | $\ell \ge 2$ |

> **On the gravitational sector.** Since the KSQC geometry is not Ricci flat, the vacuum
> perturbation equation $\delta R_{\mu\nu}=0$ cannot be imposed. The axial sector is therefore
> modelled as an **effective-source Regge–Wheeler** problem, interpreting the background as
> supported by an effective anisotropic stress tensor that carries no independent axial degree of
> freedom. The resulting spectrum should be read as that of this effective-source model, **not**
> as the full gravitational QNM spectrum of a specified quantum-gravity action.

### Two deformation-parameter normalisations

The literature uses two different normalisations of the deformation parameter, and the data files
here follow suit. They are related by

$$\mathfrak{a} = \frac{a}{M}, \qquad
a_K = \frac{a}{r_h} = \frac{\mathfrak{a}}{\sqrt{\mathfrak{a}^2 + 4}}, \qquad
\mathfrak{a} = \frac{2\,a_K}{\sqrt{1 - a_K^2}}.$$

- $\mathfrak{a} = a/M$ — **mass-normalised**; ranges over $[0,\infty)$.
- $a_K = a/r_h$ — **horizon-normalised** (used in the tables of Konoplya 2020); ranges over
  $[0,1)$. The near-extremal limit is $a_K \to 1^-$, equivalently $\mathfrak{a}\to\infty$.

A few reference values:

| $a_K$ | 0 | 0.3 | 0.5 | 0.9 | 0.99 |
|---|---|---|---|---|---|
| $\mathfrak{a}$ | 0 | 0.629 | 1.155 | 4.129 | 14.04 |

---

## 2. Method in brief

1. **Factor out the asymptotics.** The ingoing-at-the-horizon and outgoing-at-infinity QNM
   boundary conditions are removed analytically, leaving a residual radial function that is
   regular on the compact interval $[-1,1]$.
2. **Chebyshev collocation.** The regular part is expanded in Chebyshev polynomials of the first
   kind and collocated on the Chebyshev–Gauss roots grid. The residual is required to vanish at
   the collocation points.
3. **Quadratic matrix pencil.** This yields

   $$\big(M_0 + i\,\Omega\,M_1 + \Omega^2 M_2\big)\,\mathbf{a} = 0,$$

   with real $N\times N$ matrices $M_0, M_1, M_2$. The same quadratic structure holds for all
   four sectors; only the potential-dependent coefficient changes between them.
4. **Linearise and solve.** The pencil is linearised to a generalized eigenvalue problem of twice
   the dimension and solved at three spectral resolutions $N = 190, 195, 200$.
5. **Stability filter.** A candidate frequency is retained only if it can be matched consistently
   across all three resolutions (relative/absolute tolerance $10^{-6}$); the reported value is the
   average $\Omega_{\mathrm{av}} = (\Omega_{190} + \Omega_{195} + \Omega_{200})/3$, and the
   pairwise differences are stored as convergence diagnostics. This discriminates physical modes
   from spurious numerical roots.
6. **Instability diagnostic.** The *unfiltered* spectra are separately scanned for roots with
   $\mathrm{Im}\,\Omega > 0$; these are subjected to the same matching criterion, providing a
   direct check that the damped tables hide no stable growing mode.

The matrices $M_0, M_1, M_2$ are assembled symbolically in **Maple** with **200 decimal digits**;
the eigenvalue problem and all post-processing are performed in **Julia** with **220 decimal
digits** of working precision.

**Selected findings.** After accounting for the differing normalisations, the spectral results
reproduce the available time-domain, WKB, Mashhoon, and AIM data in their common ranges of
validity, with the best agreement against time-domain and AIM. The method additionally resolves
higher overtones and **purely imaginary, overdamped** branches whose spacing is close to the
surface-gravity scale $M\kappa = 1/4$ of the Schwarzschild black hole. In the near-extremal
regime the detected purely imaginary ladders approach $\Omega_n \simeq -\tfrac{i}{4}(n+1)$ in the
scalar, electromagnetic, and axial sectors, suggesting an approximate near-extremal degeneracy,
while the oscillatory spectra remain sector-dependent.

---

## 3. Repository layout

```
KS-quantum/
├── pdf/
│   └── QNM_Kazakov_Solodukhin_quantum_corrected_black_hole.pdf   # the manuscript
│
├── maple/                                  # symbolic matrix assembly & derivations (Maple)
│   ├── Matrix_Assembler_all_s_0_s_1_KS-Final Version.mw   # assembles M0/M1/M2 (scalar & EM)
│   ├── 01-CASE-s0_s1-SUPPLEMENTAL MATERIAL.mw             # scalar + electromagnetic
│   ├── 01-CASE-s0conf-SUPPLEMENTAL MATERIAL.mw            # conformally (non-min.) coupled scalar
│   └── 01-CASE-s2-SUPPLEMENTAL MATERIAL.mw                # axial gravitational sector
│
├── julia/
│   └── KS_spectral_solver_with_unstable_check.jl          # eigensolver + stability/instability filter
│
├── reports/                                # tabulated, human-readable stable-QNM tables
│   ├── s0/             (48 files)  — minimally coupled scalar
│   ├── s1/             (59 files)  — electromagnetic
│   ├── s2-axial/       (24 files)  — axial gravitational (effective source)
│   └── s0-conformal/   (10 files)  — non-minimally (conformally) coupled scalar
│
├── 200digit200l0s0/                        # convergence study — scalar monopole (s=0, ℓ=0)
│   ├── data-0/   data-5/   data-10/        #   deformation 𝔞 = 0, 5, 10
│   │   ├── eigs-190-*_full.txt             #   raw spectra at N = 190, 195, 200
│   │   ├── eigs-195-*_full.txt
│   │   ├── eigs-200-*_full.txt
│   │   ├── stable_modes.txt                #   modes surviving the 3-resolution filter
│   │   └── *.png                           #   spectrum plots (all / stable only)
│   └── …                                   #   200-digit precision, N = 200 basis
│
├── 200Digit400l0s0/                        # same study at a larger N = 400 basis (resolution check)
│   └── data-0/  data-5/  data-10/  →  qnm_plot.png, qnm_plot_stable.png
│
├── 400digit400L0S0/                        # 400-digit precision, N = 400 — parameter-sweep summary
│   ├── real_vs_a.png                       #   Re Ω vs deformation
│   └── imag_vs_a.png                       #   Im Ω vs deformation
│
├── LICENSE                                 # GNU LGPL v2.1
└── README.md
```

> **Note on `reports/s0-conformal/`.** This directory holds the **non-minimally (conformally)
> coupled** scalar data (the sector with the extra $R/6$ curvature term, assembled by
> `maple/01-CASE-s0conf-SUPPLEMENTAL MATERIAL.mw`). It was previously named `s0-min_coupled`; the
> name has been corrected so that it reflects the physics. As a cross-check, at $a_K = 0.9$ the
> fundamental mode here is $0.063952 - 0.084053\,i$ (the paper's $\Omega_c$), whereas the genuinely
> minimally coupled scalar in `reports/s0/` gives $0.051176 - 0.079180\,i$. The inner file names
> keep the `s0_` prefix because the spin is genuinely $0$.

---

## 4. File formats

### `reports/<sector>/*.txt` — tabulated stable QNMs

Each file lists the QNM frequencies for one (sector, multipole $\ell$, deformation) triple that
survive the three-resolution stability filter. A short header documents the sector and the
pencil, followed by a fixed-width table:

```
# Stable KS QNM modes
# Sector: s = 0, field = scalar, L = 0, epsilon = 1, a-tag = 0
# File tag: s0_L0_a0
# Pencil: M0 + i*Omega*M1 + Omega^2*M2 = 0
# Re/Im are printed with 16 decimal places.
# Errors are printed in scientific notation with 6 digits.
...
    k                 Re(avg)                 Im(avg)             err12             err23             err13            maxerr
-----------------------------------------------------------------------------------------------------------------------------
    5     -0.1104549390804197     -0.1048957170868810      2.539232e-18      1.285081e-18      1.761431e-18      2.539232e-18
    6      0.1104549390804197     -0.1048957170868810      2.539232e-18      1.285081e-18      1.761431e-18      2.539232e-18
   11      0.0000000000000000    -18.1166325182613148      9.089133e-05      7.632793e-04      6.723879e-04      7.632793e-04
```

| Column | Meaning |
|---|---|
| `k` | running index of the retained mode |
| `Re(avg)`, `Im(avg)` | real / imaginary part of $\Omega_{\mathrm{av}}$ (16 decimal places) |
| `err12`, `err23`, `err13` | pairwise distances between the $N=190,195,200$ estimates |
| `maxerr` | the largest of the three pairwise distances (a convergence indicator) |

Oscillatory modes appear in $\pm\,\mathrm{Re}\,\Omega$ pairs; purely imaginary overdamped modes
have `Re(avg) = 0`. A companion `*_full_precision.txt` file (regenerated by the solver, not all
checked in) carries the same data at full working precision.

**Filename convention:** `s<spin>_L<ℓ>_a<tag>.txt`, where the decimal point is written as an
underscore. The deformation `<tag>` uses one of two normalisations:

- `aK0_9` → horizon-normalised $a_K = 0.9$;
- `a1_25` → mass-normalised $\mathfrak{a} = 1.25$;
- `a0` → $\mathfrak{a} = a_K = 0$ (the Schwarzschild limit).

### Convergence-study directories (`200digit200l0s0/`, `200Digit400l0s0/`, `400digit400L0S0/`)

These document the robustness of the scalar-monopole spectrum — in particular the overdamped
purely imaginary ladder — under increasing precision (200 → 400 digits) and basis size
(N = 200 → 400). Subfolders `data-0/`, `data-5/`, `data-10/` correspond to deformation
$\mathfrak{a} = 0, 5, 10$ (Schwarzschild through near-extremal).

- `eigs-<N>-<𝔞>_full.txt` — the **raw** spectrum at resolution `N`: two whitespace-separated
  columns, `Re(Ω)` and `Im(Ω)`, at full precision (before any filtering).
- `stable_modes.txt` — the modes surviving the three-resolution match, with their matched
  partners and pairwise distances.
- `*.png` — spectrum plots in the complex $\Omega$ plane (all candidates vs. stable-only), and
  in `400digit400L0S0/`, the parameter sweeps `real_vs_a.png` / `imag_vs_a.png`.

---

## 5. Reproducing the results

The pipeline has two stages: **Maple** assembles the spectral matrices, then **Julia** solves and
post-processes them.

### Stage 1 — assemble the matrices (Maple)

Open the assembler worksheet for the desired sector:

- scalar / electromagnetic → `maple/Matrix_Assembler_all_s_0_s_1_KS-Final Version.mw`
- non-minimally coupled scalar → `maple/01-CASE-s0conf-SUPPLEMENTAL MATERIAL.mw`
- axial gravitational → `maple/01-CASE-s2-SUPPLEMENTAL MATERIAL.mw`

Set the spin `s`, the multipole `L`, and the deformation, then export the three matrices for each
resolution $N \in \{190, 195, 200\}$ at 200-digit precision. The assembler writes ASCII files
named

```
M0-<N>-s<S>_L<L>_a<TAG>.mat
M1-<N>-s<S>_L<L>_a<TAG>.mat
M2-<N>-s<S>_L<L>_a<TAG>.mat
```

where `<TAG>` is the **integer** $\mathfrak{a}\times 1000$ (e.g. `629` $\Rightarrow \mathfrak{a}=0.629$,
i.e. $a_K = 0.3$; `0` $\Rightarrow$ Schwarzschild).

### Stage 2 — solve and post-process (Julia)

The solver
[`julia/KS_spectral_solver_with_unstable_check.jl`](julia/KS_spectral_solver_with_unstable_check.jl)
reads the matrices, solves the linearised generalized eigenvalue problem in `BigFloat` (220
digits) using a generic Schur factorisation, applies the stability and instability filters, and
writes the spectra, stable-mode reports, and plots.

Install the required packages once:

```julia
import Pkg
Pkg.add(["GenericSchur", "Plots", "PlotlyJS"])   # DelimitedFiles, Printf, LinearAlgebra are stdlib
```

Choose the sector/multipole/deformation by editing the three constants near the top of the file
so that they match the Maple file tag:

```julia
const SPIN = 1               # 0 = scalar, 1 = electromagnetic, 2 = gravitational
const ELL  = 1               # multipole ℓ  (s=0: ℓ≥0, s=1: ℓ≥1, s=2: ℓ≥2)
const DEFORMATION_TAG = "629"  # integer 𝔞×1000, matching the .mat file name
```

Point the solver at the matrices and an output directory via environment variables (this avoids
editing the hard-coded default path), then run:

```bash
export KS_DATA_DIR=/path/to/maple/matrices
export KS_OUT_DIR=/path/to/output
julia julia/KS_spectral_solver_with_unstable_check.jl
```

For each run the solver writes, into `KS_OUT_DIR`:

- `spectrum_finite-N*-…​.dat` and `spectrum_candidates-N*-…​.dat` — raw and filtered spectra;
- `stable_modes-…​.txt` (+ `_details`, `_full_precision`) — the retained physical modes;
- `stable_modes_oscillatory-…​.dat` and `stable_modes_pure_imag-…​.dat` — the two branches;
- `unstable_modes-…​.txt` — the instability diagnostic (states *"No unstable mode detected."* when
  the damped tables are complete);
- `qnm_plot*-…​.png` / `.html` — complex-plane spectrum plots.

---

## 6. How to cite

If you use this code or data, please cite the accompanying paper:

> D. Batic, D. Dutykh, and M. E. Sukaiti, *Quasinormal modes of the Kazakov–Solodukhin
> quantum-corrected black hole: a spectral analysis* (2026).
> *Journal / arXiv / DOI: to be added upon publication.*

---

## 7. License

This repository is released under the **GNU Lesser General Public License, version 2.1**. See
[`LICENSE`](LICENSE) for the full text.

---

## 8. Authors and contact

- **Davide Batic** — `davide.batic@ku.ac.ae`
- **Denys Dutykh** — `denys.dutykh@ku.ac.ae`
- **Mark Essa Sukaiti** — `100064482@ku.ac.ae`

*Mathematics Department, Khalifa University of Science and Technology, PO Box 127788, Abu Dhabi,
United Arab Emirates.*
