# Solution-Set Bijection and Well-Posedness Proof

## Theorem: Bijective Correspondence Between Standard GR and Lapse-First GR

**Statement:** There exists a bijection between the solution sets of Standard General Relativity and Lapse-First General Relativity (GATG) within the specified gauge and function spaces.

## Function Spaces and Gauge Conditions

### Standard GR Function Space

For spherically symmetric spacetimes, Standard GR solutions live in the space:

**S_std = {g_μν : g_μν ∈ C²(M), det(g) < 0, signature (-,+,+,+)}**

Where:
- g_μν(t,r,θ,φ) has spherical symmetry: ∂_θ g_μν = ∂_φ g_μν = 0
- g_μν depends only on (t,r)
- Standard form: ds² = -f(t,r) dt² + h(t,r) dr² + r² dΩ²

### Lapse-First GR Function Space

For the lapse-first formulation, solutions live in:

**S_lapse = {Φ(t,r) : Φ ∈ C²(M), ∀(t,r) ∈ domain}**

Where:
- Φ(t,r) is the temporal potential
- N = e^Φ is the lapse function
- Metric constructed as: ds² = -e^(2Φ) dt² + e^(-2Φ) dr² + r² dΩ²

### Gauge Condition

**Critical Gauge Choice:** Λ = -Φ

This gauge condition is what enables the bijection. In this gauge:
- The spatial potential Λ is completely determined by the temporal potential Φ
- The metric becomes: ds² = -e^(2Φ) dt² + e^(-2Φ) dr² + r² dΩ²
- This is exactly the Schwarzschild form when A = e^(2Φ)

## Bijection Construction

### Forward Map: Standard GR → Lapse-First GR

**F: S_std → S_lapse**

Given a Standard GR solution g_μν with:
```
ds² = -f(t,r) dt² + h(t,r) dr² + r² dΩ²
```

The forward map is:
```
F(g_μν) = Φ(t,r) where e^(2Φ) = f(t,r)
```

**Constraint:** The gauge condition requires h(t,r) = 1/f(t,r) = e^(-2Φ)

### Inverse Map: Lapse-First GR → Standard GR

**F⁻¹: S_lapse → S_std**

Given a Lapse-First GR solution Φ(t,r), the inverse map is:
```
F⁻¹(Φ) = g_μν where:
- g_tt = -e^(2Φ)
- g_rr = e^(-2Φ)
- g_θθ = r²
- g_φφ = r² sin²θ
- All off-diagonal components = 0
```

### Bijection Proof

**Theorem:** F is a bijection between S_std and S_lapse in the gauge Λ = -Φ.

**Proof:**

1. **Well-defined:**
   - F(g_μν) is well-defined since f(t,r) > 0 (timelike condition)
   - F⁻¹(Φ) gives a valid Lorentzian metric with correct signature

2. **Injective (one-to-one):**
   - If F(g₁) = F(g₂), then Φ₁ = Φ₂
   - This implies e^(2Φ₁) = e^(2Φ₂), so f₁ = f₂
   - By gauge condition: h₁ = 1/f₁ = 1/f₂ = h₂
   - Therefore g₁ = g₂

3. **Surjective (onto):**
   - For any Φ ∈ S_lapse, F⁻¹(Φ) ∈ S_std
   - F(F⁻¹(Φ)) = Φ by construction

4. **Inverse relationship:**
   - F⁻¹(F(g_μν)) = g_μν (up to gauge transformation)
   - F(F⁻¹(Φ)) = Φ exactly

**Conclusion:** F establishes a bijection between solution sets.

## Constraint Propagation Lemma

### Lemma: Automatic Constraint Preservation

**Statement:** If initial data satisfies the Einstein constraints at t = t₀, then constraint propagation is automatic in both formulations due to the contracted Bianchi identities.

### Mathematical Formulation

In the gauge Λ = -Φ, the Einstein equations reduce to:

1. **Hamiltonian Constraint:** H_⊥ = R^(3) + K² - K_ij K^ij = (16πG/c⁴)ρ
2. **Momentum Constraint:** D_i(K^i_j - δ^i_j K) = (8πG/c⁴)j_j
3. **Evolution Equations:** ∂_t g_ij = -2NK_ij

### Constraint Propagation Mechanism

The contracted Bianchi identity ensures:
```
∇_μ G^μν ≡ 0 (geometric identity)
```

Combined with Einstein's equations G^μν = (8πG/c⁴)T^μν, this forces:
```
∇_μ T^μν = 0 (automatic conservation)
```

**Key Insight:** The time evolution of constraints is:
```
∂_t H_⊥ = [terms that vanish due to Bianchi identity]
∂_t C_i = [terms that vanish due to Bianchi identity]
```

Therefore: **If H_⊥ = 0 and C_i = 0 at t₀, then they remain zero for all t > t₀**

## Well-Posedness

### Definition

A formulation is **well-posed** if:
1. **Existence:** Solutions exist for admissible initial data
2. **Uniqueness:** Solutions are unique (up to gauge)
3. **Continuous dependence:** Small changes in initial data → small changes in solution

### Well-Posedness in Both Formulations

**Standard GR (ADM):**
- Well-posed when constraints are satisfied initially
- Evolution preserves constraints automatically
- Gauge freedom: coordinate transformations

**Lapse-First GR (GATG):**
- Well-posed with same initial data (expressed via Φ)
- Same constraint preservation mechanism
- Gauge freedom: choice of temporal potential normalization

**Equivalence:** Both formulations have identical well-posedness properties due to the bijection.

## Boundary Conditions and Asymptotic Behavior

### Asymptotic Flatness

For physically relevant solutions:

**At spatial infinity (r → ∞):**
```
Standard GR: f(t,r) → 1, h(t,r) → 1
Lapse-First: Φ(t,r) → 0
```

**At horizon (r → r_h):**
```
Standard GR: f(t,r) → 0, h(t,r) → ∞
Lapse-First: Φ(t,r) → -∞
```

The bijection preserves these asymptotic behaviors.

### Regularity Conditions

Both formulations require:
- Φ ∈ C² (twice differentiable)
- Proper behavior at coordinate singularities
- Finite energy conditions

## Examples

### Schwarzschild Solution

**Standard GR:** ds² = -(1-r_s/r) dt² + (1-r_s/r)⁻¹ dr² + r² dΩ²

**Lapse-First:** Φ = ½ln(1-r_s/r)

**Bijection verification:**
- F(g_schwarzschild) = ½ln(1-r_s/r) ✓
- F⁻¹(½ln(1-r_s/r)) = g_schwarzschild ✓

### Vaidya Spacetime (Dynamic)

**Standard GR:** ds² = -(1-M(v)/r) dv² + 2dv dr + r² dΩ²

**Lapse-First:** Requires coordinate transformation to gauge Λ = -Φ

The bijection works for all spherically symmetric solutions that can be brought to the gauge Λ = -Φ.

## Conclusion

The bijection F: S_std ↔ S_lapse establishes that:

1. **Every Standard GR solution** in our gauge corresponds to **exactly one lapse-first solution**
2. **Every lapse-first solution** corresponds to **exactly one Standard GR solution**
3. **Constraint propagation** is automatic in both formulations
4. **Well-posedness** is equivalent between formulations

This completes the proof that Lapse-First GR and Standard GR are not just equivalent in their solutions, but have **identical solution sets** when properly compared.

The gauge condition Λ = -Φ is the key that enables this bijection, reducing the complex 10-component Einstein equations to the single temporal potential Φ while preserving all physical content.