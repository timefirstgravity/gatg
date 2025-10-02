# Quantum Origin Module: Derivation Summary

**Reference**: Adam Snyder, "The Quantum Origin of Classical Spacetime: Gravity from Quantum Redundancy"
**DOI**: https://doi.org/10.5281/zenodo.17015383

---

## Key Claim: The CTP Kernel is DERIVED, Not Assumed

The CTP (Closed Time Path) kernel `C[Φ]` that appears in the fixed-point equation:

```
Φ = T[Φ] ≡ C ω² C[Φ]
```

is **NOT an assumption**. It is **derived from standard quantum mechanics** using the Caldeira-Leggett/Keldysh open quantum systems formalism.

---

## The Derivation Chain

### 1. Starting Point: Quantum Mechanical Bath (Appendix E.1, page 27)

**Microscopic model**:
- Large-N bosonic bath `{q_a}` with Hamiltonian `H_B = Σ_a [p²_a/(2m_a) + ½m_a ω²_a q²_a]`
- Linear coupling to lapse field: `H_int = -Σ_a g_a q_a Φ(x_a, t)`
- Bath spectral density: `J(ω) = Σ_a (g²_a/2m_a ω_a) δ(ω - ω_a)` (microphysics input)
- Initial thermal state at temperature `T_B`

**This is standard quantum mechanics** - no gravity assumptions yet.

---

### 2. Integrate Out the Bath (Appendix D.3, page 24)

**Method**: Closed-Time-Path (CTP) / Keldysh formalism

**Result**: Influence functional `S_IF[Φ⁺, Φ⁻]` (Eq. 87, page 24)

```
I[Φ⁺, Φ⁻] = exp{(i/ℏ) S_IF[Φ⁺, Φ⁻]}
          ≈ exp{(i/ℏ) ΔS[Φ] - ½ ∫∫ d⁴x d⁴x' ΔΦ(x) N(x,x') ΔΦ(x')}
```

where:
- `ΔΦ ≡ Φ⁺ - Φ⁻` (difference between forward and backward paths)
- `N(x,x')` is the **Keldysh (noise) kernel** - derived from bath integration
- `ΔS` is the dissipative (real) part

**Key Point**: `N` is **computed**, not assumed. It depends on `J(ω)` and `T_B`.

---

### 3. Generalized Langevin Equation (Eq. 99, page 28)

The influence functional gives a **generalized Langevin equation** for lapse fluctuations:

```
(-∇² + m²) δΦ(t,x) + ∫ d³x' ∫_{-∞}^t dt' Γ(x-x', t-t') ∂_{t'} δΦ(t',x') = η(t,x)
```

where:
- `Γ` is the **memory (damping) kernel** (retarded, causal)
- `η(t,x)` is the **stochastic noise** (Gaussian, zero mean)
- Both are **derived from the bath** - not free parameters

---

### 4. Fluctuation-Dissipation Theorem (Eq. 101, page 28)

The noise spectrum is **completely determined** by the damping kernel via FDT:

```
S_η(ω, k) = 2 ν(ω, T_B) Re Γ(ω, k)
```

where:
```
ν(ω, T_B) = {
    k_B T_B,                    (classical limit, high T)
    (ℏω/2) coth(ℏω/2k_B T_B),  (quantum, full FDT)
}
```

**No free parameters**: `S_η` is fixed by `(Γ, T_B)` from the bath.

---

### 5. Lapse Spectrum from Response Function (Eq. 100, page 28)

The **response function** for lapse fluctuations:

```
G_R(ω, k) = [k² + m² - iω Γ(ω, k)]^{-1}
```

The **lapse power spectrum** follows from linear response:

```
S_Φ(ω, k) = |G_R(ω, k)|² S_η(ω, k)
            = |G_R(ω, k)|² × 2 ν(ω, T_B) Re Γ(ω, k)
```

**This is the CTP noise kernel** `G^K_Φ` from Assumption A2 (page 2).

---

### 6. Capacity from Spectral Integration (Eq. 102-106, page 28)

The **capacity functional** `Ξ(x)` is the integrated lapse correlator:

```
Ξ(x) = ∫_0^T ∫_0^T C_Φ(x; t,t') dt dt'
     = ∫ (dω/2π) |W_T(ω)|² S_Φ(ω)
```

where `W_T(ω)` is the interrogation window function (rectangular: `|W_T(ω)|² = 4 sin²(ωT/2)/ω²`)

**Long-time limit** (T ≫ bath correlation time):

```
Ξ ≃ T × S_Φ(0)
```

---

### 7. Ohmic Bath Example (Eq. 107, page 29)

For an **Ohmic bath** with `Γ(ω,k) ≈ η` (real, constant at low ω):

```
S_Φ^{loc}(0) = ν(0, T_B) × η / (4πm)
             = (k_B T_B / ℏ) × η / (4πm)
```

This gives:

```
Ξ ≃ T × (k_B T_B / ℏ) × η / (4πm)
```

**All parameters are physical**:
- `T`: Interrogation time (seconds)
- `T_B`: Bath temperature (Kelvin)
- `η`: Ohmic damping strength (dimensionless)
- `m`: Screening mass (from redundancy functional curvature, Eq. 49, page 14)

---

## Summary: No Circularity, No Assumptions

The derivation proceeds as:

```
Microscopic Bath (J(ω), T_B)
    ↓ [CTP integration]
Influence Functional S_IF
    ↓ [Extract kernels]
Langevin Equation: Γ(ω), η(t)
    ↓ [FDT]
Noise Spectrum: S_η = 2ν Re Γ
    ↓ [Linear Response]
Lapse Spectrum: S_Φ = |G_R|² S_η
    ↓ [Time integration]
Capacity: Ξ = ∫∫ C_Φ dt dt'
    ↓ [Fixed-point equation]
Emergent Lapse: Φ = C ω² Ξ
```

**Every step is standard quantum open systems theory**:
- Caldeira-Leggett (1983) for quantum dissipation
- Keldysh (1965) for non-equilibrium CTP formalism
- Callen-Welton (1951) / Kubo (1966) for FDT

**No free parameters beyond microphysics**: The bath spectral density `J(ω)` and temperature `T_B` fix everything.

---

## Where the Code Implements This

### `ctp_kernel.py`
- Builds the kernel operators C (local, quasi-local, screened)
- Implements Yukawa Green function `G(r) = e^{-mr}/(4πr)` for screened case
- Computes spectral properties

### `capacity_functional.py`
- `compute_capacity_from_kernel()`: Full physics chain (Steps 1-6 above)
- `compute_ctp_noise_spectrum()`: FDT implementation (Step 4)
- `verify_fdt_relation()`: Check S_Φ = |G_R|² S_η

### `fixed_point.py`
- Constructs operator `T[Φ] = C ω² C[Φ]`
- Verifies Banach contraction: `||T|| < 1`
- Iterates to fixed point

### `physical_parameters.py`
- Benchmark parameters from paper (optical clock, Earth geoid)
- Solves for coupling that satisfies contraction constraint

---

## References to Paper

| Code Module | Paper Section | Equation | Page |
|------------|---------------|----------|------|
| CTP kernel derivation | Appendix E.1-E.2 | 87, 99-101 | 24, 27-28 |
| Capacity Ξ from spectrum | Appendix E.3 | 102-106 | 28 |
| Ohmic bath result | Appendix E.3 | 107 | 29 |
| FDT relation | Section 5.1 | 52-53 | 14 |
| Fixed-point equation | Section 3 | 15, 38 | 6, 11 |
| Screening length | Section 5.1 | 49 | 14 |
| Benchmark predictions | Section 8 | — | 17-18 |

---

## The Bottom Line

**The CTP kernel is derived from standard quantum mechanics using well-established methods (Caldeira-Leggett, Keldysh, FDT).**

**It is NOT an assumption, NOT circular, and has NO free parameters beyond the physical bath properties `(J(ω), T_B, m)`.**

**This module implements that derivation with full computational rigor.**
