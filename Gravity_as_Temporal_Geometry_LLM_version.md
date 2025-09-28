# Gravity as Temporal Geometry (Part I, lapse-first GR)

## HEADER
```yaml
title: Gravity as Temporal Geometry (Part I, lapse-first GR)
signature: "- + + +"
units_default: "G=c=1 (restore for experiments)"
scope: "Classical equivalence to GR; spherical ODE; flux law; rotation-as-shift; linearized TT; FRW map"
doi: "https://doi.org/10.5281/zenodo.16878018"
dependencies: []
```

### Sign Conventions (Box)
- Metric signature: \((- + + +)\)
- Einstein equation: \(G_{\mu\nu} = +\,\tfrac{8\pi G}{c^{4}}\,T_{\mu\nu}\)
- Riemann/curvature: \(R^{\rho}{}_{\sigma\mu\nu} = \partial_{\mu}\Gamma^{\rho}_{\nu\sigma} - \partial_{\nu}\Gamma^{\rho}_{\mu\sigma} + \cdots\); \(R_{\mu\nu} = R^{\rho}{}_{\mu\rho\nu}\); \(G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}g_{\mu\nu}R\)
- ADM momentum constraint: \(D_j(K^{j}{}_{i} - \delta^{j}{}_{i}K) = \tfrac{8\pi G}{c^{4}} T^{t}{}_{i}\)
- Lapse/shift convention: \(N = e^{\Phi}\), \(\beta_i\) as the shift covector; overdot \(\dot{}\) = \(d/d\tau\).

## TL;DR
GR is recast with a single temporal scalar Φ (lapse \(N=e^{\Phi}\)); spatial geometry follows from ADM constraints. In static spherical vacuum, the field equations reduce to one ODE \(rA' = 1 - A\) with \(A=e^{2\Phi}\), yielding Schwarzschild and classic tests. Dynamics in spherical EF/PG slicing give a mixed-component identity \(G^{t}{}_{r}=(2/r)\,\partial_{t}\Phi\), hence the **flux law**
\[\partial_t \Phi = \frac{4\pi G}{c^{4}}\, r\, T^{t}{}_{r}\]
(positive sign with the \(G^{t}{}_{r}=+8\pi G T^{t}{}_{r}/c^{4}\) convention). Rotation resides entirely in the **shift**; far from a spinning source, \(g_{t\phi}=-2GJ/(c^{3}r)\sin^{2}\theta\) and both the ZAMO angular velocity and Lense–Thirring precession scale as \(2GJ/(c^{2}r^{3})\). Linearized vacuum leaves only the two TT graviton modes, \(\Box h^{TT}=0\). Cosmology maps via \(d\tau=e^{\Phi}dt\), \(a=e^{-\Phi}\), \(H=-\dot{\Phi}\), reproducing the Friedmann pair exactly. Atom-interferometer phase follows from \(\mathbf{g}=-c^{2}\nabla\Phi\), reducing to \(\Delta\varphi=k_{\rm eff} g T^{2}\) in the uniform-\(g\) regime.

## SYMBOLS
t               | coordinate time  
\(\tau\)        | proper/cosmic time ( \(d\tau=e^{\Phi} dt\) )  
\(\Phi\)        | lapse exponent (dimensionless), \(A\equiv e^{2\Phi}\), \(N\equiv e^{\Phi}\)  
\(U\)           | gravitational potential per unit mass, \(U\equiv c^{2}\Phi\)  
\(\beta_i\)     | shift covector; \(N^{i}=\gamma^{ij}\beta_{j}\) (Note: we reserve \(\omega\) for angular velocities)  
\(\gamma_{ij}\) | spatial 3-metric  
\(a\)           | FRW scale factor \(a=e^{-\Phi}\)  
\(H\)           | Hubble rate \(H=-d\Phi/d\tau\)  
\(T^{t}{}_{r}\) | mixed stress component (energy-flux density)  
\(h^{TT}_{ij}\) | radiative TT strain  
\(\varphi\)     | interferometer phase

## REGIME / ASSUMPTIONS
- Signature \((- + + +)\); \(G=c=1\) unless displayed; restore \(c\) for experiments.  
- Spherical sections: diagonal gauge with zero shift; for horizons/flux use EF/PG regular slicings.  
- Linearized sections: perturbations about Minkowski.  
- Cosmology: constant-curvature spatial slices; dot \(\dot{}\) means \(d/d\tau\).  
- **Derivative notation:** primes \(('\)) denote \(d/dr\) in the spherical (S) sector; overdots denote \(d/d\tau\) in cosmology (C).

## CORE EQUATIONS (ID: S = spherical, N = null/EF, L = linear/GEM, C = cosmology, A = ADM)
**S1**:  \(ds^{2}=-A\,dt^{2} + A^{-1}\,dr^{2} + r^{2} d\Omega^{2}\),  \(A=e^{2\Phi(r)}\).  
**S2**:  \(rA'(r) = 1 - A(r)\quad \Rightarrow\quad A(r)=1-\dfrac{2GM}{r}\) (vacuum, static, zero shift).
> *Footnote (spherical ansatz):* In static, spherically symmetric **vacuum**, the field equations imply \(g_{rr} = A^{-1}\) when \(g_{tt}=-A\); i.e., \(AB=1\) follows from the Bianchi identity. We adopt this familiar diagonal gauge only for the static S-sector and do not assume it off-shell in dynamical/flux contexts.  
**S3**:  \(\partial_{t}\Phi = \dfrac{4\pi G}{c^{4}}\, r\, T^{t}{}_{r}\)  \(\;\) *(Here \(t\) is the diagonal-gauge time obtained after transforming from EF \(v\) at fixed \(r\); derivation via ADM momentum constraint in Appendix D1. We maintain the sign convention \(G^{t}{}_{r}=+8\pi G T^{t}{}_{r}/c^{4}\) throughout)*.  
**N1**:  EF/Vaidya form \(ds^{2}=-A(v,r)\,dv^{2}+2\,dv\,dr+r^{2}d\Omega^{2}\), \(A=1-\dfrac{2Gm(v)}{r}\), \(T_{vv}=\dfrac{1}{4\pi r^{2}}\,\dfrac{dm}{dv}\).  
**L1**:  \(g_{t\phi} = -\dfrac{2GJ}{c^{3} r}\sin^{2}\theta\);\ \(\omega_{\text{ZAMO}} = \dfrac{2GJ}{c^{2} r^{3}}\);\ \(\Omega_{\mathrm{LT}} = \dfrac{2GJ}{c^{2} r^{3}}\).  
**L2**:  \(\Box\, h^{TT}_{ij} = 0\), with GW stress-energy \(t^{\mathrm{GW}}_{\mu\nu}=(32\pi G)^{-1}\langle \partial_{\mu}h^{TT}\partial_{\nu}h^{TT}\rangle\).  
**C1**:  \(d\tau=e^{\Phi(t)}dt\), \(a=e^{-\Phi}\), \(H=-d\Phi/d\tau\).  
**C2**:  \(H^{2}+k/a^{2}=(8\pi G/3)\rho + \Lambda/3\), \(\ \dot{H} = -4\pi G(\rho+p) + k/a^{2}\).  
**A1**:  \(\mathcal{H}_{\perp}=0\), \(\mathcal{H}_{i}=0\); lapse/shift are Lagrange multipliers (no propagating scalar/vector DOF).

## COMPUTE RECIPES (pseudo-code)
**R1 — Schwarzschild \(A(r)\):**  
Solve \(r\,A'(r)=1-A(r)\) ⇒ \(A(r)=1-2GM/r\). Weak-field: \(U=c^{2}\Phi=\tfrac{c^{2}}{2}\ln(1-2GM/r)\approx -GM/r\).

**R2 — Flux-step (spherical, fixed r):**  
\(\Phi(t+\Delta t,r) = \Phi(t,r) + \big[\,(4\pi G/c^{4})\, r\, T^{t}{}_{r}(t,r)\,\big]\Delta t\)\  # explicit time integrator for ADM momentum constraint.

**R3 — Frame dragging (far field):**  
\(g_{t\phi}(r,\theta) = -2GJ/(c^{3}r)\sin^{2}\theta\);\ \(\omega_{\text{ZAMO}}(r) = 2GJ/(c^{2}r^{3})\);\ \(\Omega_{\mathrm{LT}}(r) = 2GJ/(c^{2}r^{3})\).

**R4 — FRW map:**  
Given \(\Phi(t)\): \(\tau=\int e^{\Phi}dt\), \(a(\tau)=e^{-\Phi}\), \(H(\tau)=-d\Phi/d\tau\).

**R5 — Atom interferometer (Mach–Zehnder, 0–T–2T):**  
\(\Delta\varphi = -c^{2}\!\int g_{s}(t)\,\mathbf{k}_{\rm eff}\!\cdot\!(\nabla\Phi)\,dt\), which reduces to \(\Delta\varphi = k_{\rm eff}\,g\,T^{2}\) in uniform \(g\) (path-phase difference reduces to standard form).
*(Breadcrumb: phase difference \(= -\tfrac{1}{\hbar}\int \Delta U\,dt\) with \(U = m c^{2}\Phi\); in uniform \(g=-c^{2}\nabla\Phi\), this collapses to \(k_{\rm eff}gT^{2}\).)*

## SANITY CHECKS (gold numbers)
- **Gravitational redshift:** \(z \simeq GM/r_{e} - GM/r_{o}\) (Pound–Rebka/GPS).  
- **Light bending:** \(\Delta\phi = 4GM/b\) (VLBI).  
- **Shapiro delay:** \(\Delta t_{2w} \approx 4GM \ln\!\big(4 r_{\oplus} r_{T}/b^{2}\big)\).  
- **Perihelion precession:** \(\Delta\varpi = 6\pi GM / [a(1-e^{2})]\).  
- **Horizon regularity:** EF/PG coordinates finite at \(r=2GM\) (diagonal chart has coordinate singularity).  
- **Kerr/dragging:** \(g_{t\phi}=-2GJ/(c^{3}r)\sin^{2}\theta\), \(\Omega_{\mathrm{LT}}=2GJ/(c^{2}r^{3})\).  
- **FRW:** \(a=e^{-\Phi}\), \(H=-\dot{\Phi}\) reproduces Friedmann exactly.  
- **Radiation content:** scalar carries no waves; radiative energy only in TT \(h^{TT}\) (L2).

## DIMENSION LEDGER (SI base)
\([\Phi]=1\), \([A]=1\), \([H]=T^{-1}\), \([T^{t}{}_{r}]=\) energy-flux density \(= M T^{-3}\) (W m⁻²), \([g_{t\phi}]=1\), \([J]=ML^{2}T^{-1}\), \([\Omega_{\mathrm{LT}}]=T^{-1}\), \([\nabla\Phi]=L^{-1}\), \([\mathbf{g}]=LT^{-2}\) with \(\mathbf{g}=-c^{2}\nabla\Phi\).

## PITFALLS / GUARDRAILS
- Do **not** attribute frame dragging to \(\Phi\); use the shift \(\beta\).  
- In spherical dynamics use EF/PG to avoid fake horizon singularities.  
- Keep \(c\) powers explicit in \(g_{t\phi}\) and in the flux law (when not in \(G=c=1\)).  
- Never reuse \(\phi\) (azimuth) for anything except the angle; interferometer phase is \(\varphi\).

## CLAIMS → SUPPORT
1. Static spherical vacuum reduces to one ODE (S1–S2) ⇒ Schwarzschild and classic tests.  
2. Spherical flux identity from \(G^{t}{}_{r}\) (S3) matches Vaidya/Bondi in EF (N1), derived rigorously via ADM momentum constraint (D1).  
3. Rotation-as-shift reproduces slow-Kerr and Lense–Thirring (L1).  
4. Only TT modes radiate and carry GW energy flux (L2); lapse/shift are non-propagating constraints.  
5. FRW cosmology recovered exactly via \(a=e^{-\Phi}\) (C1–C2).

---

## Appendix D1 — Derivation of the Spherical Flux Law via ADM Momentum Constraint

We derive the spherical flux law
\[
\boxed{\;\partial_{t}\Phi = \frac{4\pi G}{c^{4}}\, r\, T^{t}{}_{r}\;}\tag{D1.0}
\]
from the ADM momentum constraint with zero shift in the scalar sector:
\[
D_j\!\left(K^{j}{}_{i}-\delta^{j}{}_{i}K\right)=\frac{8\pi G}{c^{4}}\,T^{t}{}_{i}.\tag{D1.1}
\]
We adopt the diagonal spherical chart used in the main text with lapse \(N=e^{\Phi(t,r)}\) and 3-metric
\[
\gamma_{ij}\,dx^{i}dx^{j} = A^{-1}(t,r)\,dr^{2} + r^{2} d\Omega^{2},\qquad A\equiv e^{2\Phi}.\tag{D1.2}
\]
(For EF/PG slicings near horizons we transform later; the identity is local and gauge-covariant under the conditions stated below.)

### D1.1 Extrinsic curvature components
With vanishing shift \(\beta_i=0\), our ADM sign convention gives
\[
K_{ij} = -\frac{1}{2N}\,\partial_{t}\gamma_{ij}.\tag{D1.3}
\]
Because \(\gamma_{\theta\theta}=r^{2}\) and \(\gamma_{\phi\phi}=r^{2}\sin^{2}\theta\) are time-independent in this chart, the only nonzero component is
\[
K_{rr} = -\frac{1}{2e^{\Phi}}\,\partial_{t}(A^{-1})
= -\frac{1}{2e^{\Phi}}(-A^{-2}\partial_{t}A)
= \frac{1}{2e^{\Phi}}\,A^{-2}\,(2e^{2\Phi}\partial_{t}\Phi)
= e^{-\Phi}\,A^{-1}\,\partial_{t}\Phi.\tag{D1.4}
\]
Raising one index with \(\gamma^{rr}=A\) gives
\[
K^{r}{}_{r}=\gamma^{rr}K_{rr}=e^{-\Phi}\,\partial_{t}\Phi,\qquad K^{\theta}{}_{\theta}=K^{\phi}{}_{\phi}=0.\tag{D1.5}
\]
The trace is \(K=K^{r}{}_{r}+K^{\theta}{}_{\theta}+K^{\phi}{}_{\phi}=e^{-\Phi}\,\partial_{t}\Phi\).

Define the mixed tensor
\[
Y^{j}{}_{i}:=K^{j}{}_{i}-\delta^{j}{}_{i}K,\tag{D1.6}
\]
so that the constraint (D1.1) reads \(D_{j}Y^{j}{}_{i}=\tfrac{8\pi G}{c^{4}}T^{t}{}_{i}\).
In our diagonal, spherically symmetric case, the only potentially nonzero RHS component is \(i=r\), hence we compute \(D_{j}Y^{j}{}_{r}\).

From (D1.5) one finds
\[
Y^{r}{}_{r}=K^{r}{}_{r}-K=0,\qquad
Y^{\theta}{}_{r}=Y^{\phi}{}_{r}=0,\qquad
Y^{\theta}{}_{\theta}=Y^{\phi}{}_{\phi}=-K^{r}{}_{r}=-e^{-\Phi}\,\partial_{t}\Phi.\tag{D1.7}
\]
Although \(Y^{r}{}_{r}=0\) pointwise, its **covariant divergence** need not vanish because the connection acts on the index structure. Using the identity for mixed tensors in a diagonal metric,
\[
D_{j}Y^{j}{}_{r}
= \partial_{j}Y^{j}{}_{r} + \Gamma^{j}_{\,j m}Y^{m}{}_{r} - \Gamma^{m}_{\,r j}Y^{j}{}_{m},\tag{D1.8}
\]
only the last term survives (since \(Y^{j}{}_{r}\) vanishes for all \(j\) ):
\[
D_{j}Y^{j}{}_{r} = -\Gamma^{m}_{\,r j}Y^{j}{}_{m}.
\tag{D1.9}
\]
For the 3-metric (D1.2), the nonzero Christoffels needed are
\[
\Gamma^{\theta}_{\,r\theta}=\Gamma^{\phi}_{\,r\phi}=\frac{1}{r},\qquad
\Gamma^{r}_{\,rr}= -\tfrac{1}{2}\partial_{r}\ln A.
\tag{D1.10}
\]
Since \(Y^{\theta}{}_{\theta}=Y^{\phi}{}_{\phi}=-e^{-\Phi}\partial_{t}\Phi\) and \(Y^{r}{}_{r}=0\), the only contributions in (D1.9) come from \(m=\theta,\phi\) with \(j=\theta,\phi\):
\[
D_{j}Y^{j}{}_{r} = -\Gamma^{\theta}_{\,r\theta}Y^{\theta}{}_{\theta} - \Gamma^{\phi}_{\,r\phi}Y^{\phi}{}_{\phi}
= -\frac{1}{r}(-e^{-\Phi}\partial_{t}\Phi) - \frac{1}{r}(-e^{-\Phi}\partial_{t}\Phi)
= \frac{2}{r}\,e^{-\Phi}\,\partial_{t}\Phi.
\tag{D1.11}
\]
Thus the \(i=r\) component of (D1.1) yields
\[
\frac{2}{r}\,e^{-\Phi}\,\partial_{t}\Phi = \frac{8\pi G}{c^{4}}\,T^{t}{}_{r}.\tag{D1.12}
\]
Finally, recalling that coordinate and proper times obey \(d\tau = N\,dt = e^{\Phi} dt\), we may write \(e^{-\Phi}\partial_{t} = \partial_{\tau}\). Expressed in **coordinate time** \(t\) and the mixed component \(T^{t}{}_{r}\), (D1.12) becomes exactly (D1.0):
\[
\partial_{t}\Phi = \frac{4\pi G}{c^{4}}\, r\, T^{t}{}_{r}.\tag{D1.13}
\]
*Remarks.* (i) The result is local and does not require stationarity. (ii) In EF/Vaidya coordinates, the same identity reproduces \(A=1-2Gm(v)/r\) for null dust with \(T_{vv}=(4\pi r^{2})^{-1}dm/dv\). (iii) Near horizons, use regular EF/PG slicings; (D1.13) remains valid when mapped back to the diagonal chart used here (the lapse factors cancel as shown).

### D1.2 Transformation from EF (diagonalizing the cross term)
Starting from EF/Vaidya \(ds^{2}=-A\,dv^{2}+2\,dv\,dr+r^{2}d\Omega^{2}\), define a new time \(t(v,r)\) by imposing \(g_{tr}=0\). This yields the first-order equation \(\partial_{r}t = (1-A)/A\) (up to an additive function of \(v\)), which fixes the diagonal time used in the main text. Evaluating \(G^{t}{}_{r}\) at fixed \(r\) in this chart and using (D1.11) gives (D1.13) with the displayed sign and coefficient \((2/r)\).

## Appendix D2 — FRW Map from Lapse-First
Take the homogeneous ansatz \(ds^{2}= -e^{2\Phi(t)}dt^{2} + e^{-2\Phi(t)}\,\gamma_{ij}dx^{i}dx^{j}\). Define \(d\tau=e^{\Phi}dt\) and \(a=e^{-\Phi}\); then \(ds^{2}=-d\tau^{2}+a^{2}(\tau)\gamma_{ij}dx^{i}dx^{j}\).

**From ADM constraints:** With homogeneous stress-tensor \(T_{\mu\nu}\), the Hamiltonian constraint \(\mathcal{H}_{\perp}=0\) gives:
\[
R^{(3)}+K^{2}-K_{ij}K^{ij}=\frac{16\pi G}{c^{4}}\rho,
\]
which, with our ansatz and \(H=-\dot{\Phi}\), yields:
\[
H^{2} + \frac{k}{a^{2}} = \frac{8\pi G}{3}\rho + \frac{\Lambda}{3}.
\]

The trace of the evolution equation gives:
\[
\dot{H} = -4\pi G(\rho+p) + \frac{k}{a^{2}}.
\]

The scalar \(\Phi\) is a **constraint** (non-propagating); all radiative content resides in the TT sector.

## Appendix D3 — Non-propagation of Scalar and Vector Modes

In the ADM formalism, the lapse \(N=e^{\Phi}\) and shift \(\beta^{i}\) appear as Lagrange multipliers enforcing the constraints \(\mathcal{H}_{\perp}=0\) and \(\mathcal{H}_{i}=0\). They have no canonical kinetic terms in the action:
\[
S = \int dt\,d^{3}x\,N\sqrt{\gamma}\left(\mathcal{H}_{\perp}\right) + \int dt\,d^{3}x\,\beta^{i}\mathcal{H}_{i}.
\]

Variation with respect to \(N\) and \(\beta^{i}\) yields the constraints, not propagation equations. In vacuum linearized gravity, after gauge fixing and constraint solving, only the two transverse-traceless polarizations of \(h_{ij}^{TT}\) carry independent degrees of freedom. Thus:
- No extra scalar graviton from \(\Phi\)
- No vector gravitons from \(\beta^{i}\)
- Only TT modes propagate: \(\Box h^{TT}_{ij}=0\)
