# Wave equation on a disc — a spectral MATLAB solver

A solver for the 2D wave equation on a disc, built on the exact Bessel
eigenbasis rather than on a mesh, plus a Phase 2 framework that carries the
same machinery to conformally mapped domains.

```
w_tt = c² Δw            in  D = { x² + y² < R² }
w    = 0                on  ∂D        (Dirichlet, default)
∂w/∂n = 0               on  ∂D        (Neumann, toggle)
w(·,0) = w₀ ,  w_t(·,0) = v₀          (smooth; Gaussian by default)
```

## Quickstart

```matlab
addpath('waveDisc')

sol = wave_solver_disc;                 % solve with the defaults
plot_solution(sol, 2.4)                 % one snapshot
sol.eval(0.1, 0.2, [0 1 2])             % w at a point, at three times

demo_wave_disc                          % all the figures and the animation
demo_wave_polygon                       % the Phase 2 figure
test_wave_solver_disc                   % 15-check verification suite
test_wave_solver_polygon                % 8-check verification of Phase 2
```

## How it works

In polar coordinates the Dirichlet eigenfunctions of −Δ on the disc are

```
φ_{m,n}(r,θ) = J_m(λ_{mn} r) · {cos mθ, sin mθ} ,   λ_{mn} = j_{m,n}/R
```

with `j_{m,n}` the n-th positive zero of `J_m` (of `J_m'` for Neumann, plus
the constant mode). Expanding `w = Σ α_k(t) φ_k` turns the PDE into K
**decoupled** harmonic oscillators `α_k'' + ω_k² α_k = 0`, `ω_k = c λ_k`.

The pipeline is therefore:

1. **`bessel_zeros`** — zeros of `J_m` and `J_m'`, by bracketing sign
   changes and polishing with `fzero`. Robust at every order; cached.
2. **`disc_eigenmodes`** — every mode with `λ ≤ LambdaMax`, with closed-form
   L² norms. The truncation is a disc in spectral space, not a box in
   `(m,n)`, which is the right balance for an isotropic problem.
3. **Projection** — Gauss–Legendre in `r` (`gauss_legendre`, Golub–Welsch,
   no toolbox) tensored with a periodic trapezoid rule in `θ` evaluated by
   FFT. Both are spectrally accurate for smooth data.
4. **Time** — closed form per mode by default; `ode45` on the same modal
   system is available through `'TimeScheme','ode45'`.
5. **Synthesis** — the field on a polar grid is recovered by one inverse
   FFT in `θ` per time step, so plotting stays cheap even with thousands of
   modes.

Nothing needs a MATLAB toolbox: only `besselj`, `fzero`, `eig`, `fft` and
`ode45`, all of which are in base MATLAB. The code also runs unmodified in
GNU Octave.

### Why the closed form is the default

The brief asked for a standard ODE solver, and `'TimeScheme','ode45'` gives
exactly that on the modal system. It is not the default because for this
particular system the closed form is both exact and far cheaper — figure 5b
shows `ode45` drifting to ~3·10⁻⁸ in energy over `t ∈ [0,20]` where the
closed form holds 10⁻¹⁶. The `ode45` path is kept because it is the hook for
anything the closed form cannot absorb: damping, forcing, or the
variable-coefficient problem of Phase 2.

### The one thing to watch: the boundary

A Gaussian never vanishes *identically* on `r = R`, so under Dirichlet
conditions the expansion of `w₀` converges only down to the level by which
the tail violates the boundary condition. The solver measures this and
warns:

```
  projection   : captures 1 - 1.6e-09 of ||w0||_L2
  BC residual  : 1.6e-09 of peak (w on dD)
```

Keep the pulse a few `σ` clear of the rim and both numbers stay negligible.
Figure 5a shows the convergence curve flattening onto exactly that floor.

## Files

| file | role |
|---|---|
| `wave_solver_disc.m` | **Phase 1 main solver** |
| `wave_solver_polygon.m` | **Phase 2 framework**, conformal mapping |
| `initial_condition.m` | Gaussian, ring, two bumps, single eigenmode, custom |
| `plot_solution.m` | one snapshot: flat, 3D relief or contour |
| `animate_wave.m` | GIF/video, rasterised directly, no `getframe` |
| `disc_eigenmodes.m` | the eigenbasis and its exact norms |
| `bessel_zeros.m`, `dbesselj.m` | zeros of `J_m`, `J_m'` |
| `gauss_legendre.m` | quadrature nodes and weights |
| `wave_colormap.m`, `wave_colorbar.m`, `wave_colorstrip.m` | colour |
| `save_fig.m` | PNG output across releases |
| `demo_wave_disc.m`, `demo_wave_polygon.m` | the figures |
| `test_wave_solver_disc.m`, `test_wave_solver_polygon.m` | verification |

## Figures produced by `demo_wave_disc`

| file | content |
|---|---|
| `fig1_evolution.png` | eight snapshots of the front and its reflection |
| `fig2_surfaces.png` | the same field as 3D relief |
| `fig3_eigenmodes.png` | gallery of the 20 lowest Dirichlet eigenfunctions |
| `fig4_spectrum.png` | modal energy, spectral decay, Weyl's law |
| `fig5_verification.png` | convergence, conservation, exactness of the BCs |
| `fig6_spacetime.png` | space-time diagrams with the characteristic overlaid |
| `fig7_dirichlet_vs_neumann.png` | the boundary-condition toggle |
| `fig8_polygon.png` | Phase 2, from `demo_wave_polygon` |
| `wave_disc.gif` | animation |

## Verification

`test_wave_solver_disc` runs 15 checks; the ones that matter most:

| check | error |
|---|---|
| Bessel zeros against tabulated values | 5·10⁻¹⁵ |
| eigenmode norms against brute-force quadrature | 1·10⁻¹⁴ |
| single mode against its closed form | 5·10⁻¹⁵ |
| Dirichlet `w = 0` on `∂D` | 2·10⁻¹⁵ |
| energy conservation | 1·10⁻¹⁶ |
| mass conservation (Neumann) | 8·10⁻¹⁵ |
| `ode45` against the closed form | 2·10⁻⁹ |
| **`w_tt − c²Δw` by finite differences** | 9·10⁻⁶ |
| projection error, `λ_max` 12 → 36 | 0.18 → 4·10⁻⁷ |

The finite-difference residual is the one that treats the solver as a black
box: it differentiates the returned solution and puts it back into the PDE.

`test_wave_solver_polygon` runs 8 more. The two exact ones pin the whole
conformal chain: under the identity map the weighted eigenproblem must
return the Bessel spectrum of the disc (error 5·10⁻¹⁵), and under `z = s·ζ`
it must return that spectrum divided by `s` (error 3·10⁻¹⁵).

## Phase 2 status

Working, but experimental, and with one piece of mathematics that is easy to
get wrong. **A conformal map does not turn the wave equation on a polygon
into the wave equation on the disc.** Under `z = f(ζ)`,

```
Δ_z = |f'(ζ)|⁻² Δ_ζ      ⟹      ρ(ζ) w_tt = c² Δ_ζ w ,   ρ = |f'|²
```

so the pulled-back problem has a variable coefficient. Only Laplace and
Helmholtz are conformally clean. What survives is separation of variables
against a **weighted** eigenproblem `−Δφ = μ² ρ φ`, which
`wave_solver_polygon` assembles in the Bessel basis (`S v = μ² M v`, with
`S` diagonal and `M` the weighted mass matrix) and diagonalises. Because
`ρ` is the Jacobian, `∫ρ u v` is exactly the physical L²(Ω) product, so no
further change of measure appears anywhere.

Known limits, in the order they will bite:

* `M` is dense, at `O(K² · N_quad)` work, so the affordable basis is a few
  hundred modes rather than the few thousand Phase 1 handles comfortably.
  A fast transform in `θ` would be the obvious next step.
* Schwarz–Christoffel maps crowd badly for elongated polygons, and `ρ` then
  spans orders of magnitude. Reasonable aspect ratios only.
* Reentrant corners make `∇w` singular, so convergence there is algebraic
  rather than spectral.
* The Neumann path goes through the same solve and keeps the zero mode, but
  has had much less exercise than Dirichlet.

For a true polygon, put Driscoll's Schwarz–Christoffel toolbox on the path
and pass vertices:

```matlab
addpath('/path/to/sc-toolbox')
sol = wave_solver_polygon('Vertices', [0, 1, 1+0.7i, 0.3+1i, 0.9i], 'T', 4);
plot_solution(sol, 1.5)
```

Without the toolbox, any analytic conformal map of the disc works and needs
no extra dependency:

```matlab
sol = wave_solver_polygon('MapFun', @(z) z + 0.22*z.^4, ...
                          'MapDiff', @(z) 1 + 0.88*z.^3, 'LambdaMax', 26);
```

## Notes on the committed artefacts

The figures and the animation in `figures/` were generated by the scripts
here, running under GNU Octave 8.4 with the Qt/OpenGL backend. Two details
differ from a MATLAB run:

* `save_fig` uses `exportgraphics` where it exists and falls back to
  `print`; Octave takes the fallback, so MATLAB output will be slightly
  crisper.
* Octave's `imwrite` writes GIFs without LZW compression, and it has no
  `VideoWriter`. The committed `wave_disc.gif` was therefore re-encoded
  losslessly with `gifsicle -O3` after `animate_wave` produced it, and
  `wave_disc.mp4` was encoded with `ffmpeg` from the same frames. On MATLAB,
  `animate_wave` writes both directly, and no post-processing is needed.
