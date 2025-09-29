# Quantum Witness Module

Implements quantum signature detection protocols for gravitational dephasing
in clock networks, based on GATG Paper V.
This module detects quantum non-commutativity through time-reversal asymmetry
in clock dephasing experiments.

## Core Concept

The **Q1 commutator witness** detects quantum signatures by measuring the
difference in visibility between time-forward and time-reversed clock protocols:

- **Classical systems**: `V_forward = V_reverse` → `Δ(-ln V) = 0`
- **Quantum systems**: `V_forward ≠ V_reverse` → `Δ(-ln V) ≠ 0`

## Module Components

### `spectral_core.py`

Fundamental spectral analysis functions:

- `build_psd_lapse_fluctuations()` - Power spectral density S_Φ(f) from GATG flux law
- `build_antisymmetric_spectrum()` - Quantum antisymmetric spectrum S⁻_Φ(f)
- `compute_total_filter()` - Combined filter F(ω) = Y(ω) × G(ω;L)
- `verify_spectral_properties()` - Mathematical verification of spectral functions

### `filters.py`

Clock sequence and baseline filters:

- `ramsey_filter()` - Ramsey interferometry sequence Y_R(ω)
- `hahn_echo_filter()` - Hahn echo sequence Y_E(ω) with noise suppression
- `cpmg_filter()` - CPMG sequence for enhanced coherence
- `baseline_response_function()` - Baseline factor G(ω;L) for separated clocks

### `commutator_witness.py`

Core quantum detection protocol:

- `compute_q1_witness()` - Q1 commutator witness: `Δ(-ln V) = 2 Im ∫ df Ω_A Ω_B* S⁻_Φ`
- `compute_visibility_difference()` - Visibility asymmetry analysis
- `verify_classical_limit()` - Confirms classical systems give zero witness
- `compute_time_reversal_protocol()` - Time-reversal implementation

### `cumulants.py`

Higher-order statistical analysis:

- `compute_second_cumulant()` - K₂ (Gaussian dephasing)
- `compute_third_cumulant()` - K₃ (non-Gaussian signatures)
- `compute_fourth_cumulant()` - K₄ᶜ (connected fourth-order)
- `compute_snr_lines()` - Signal-to-noise ratio calculations

## Mathematical Foundation

**Q1 Witness Formula**:

```markdown
Δ(-ln V) = 2 Im ∫ df Ω_A(f) Ω_B*(f) S⁻_Φ(f)
```

Where:

- `Ω_A(f), Ω_B(f)` are clock filter functions
- `S⁻_Φ(f)` is the antisymmetric lapse spectrum (pure imaginary, odd in frequency)
- The integral isolates quantum non-commutativity signatures

**Physical Regularization**: Uses Lorentzian cutoff `f² → f² + f_c²` to prevent low-frequency divergences.

## Key Features

✓ **Rigorous symbolic computation** with SageMath
✓ **No fake data** - all computations use real mathematical results
✓ **Classical limit verification** - confirms `ħ → 0` gives zero witness
✓ **Dimensional consistency** checking throughout
✓ **Multiple clock sequences** (Ramsey, Hahn echo, CPMG)
✓ **Baseline configurations** (co-located, one-way, two-way)

## Usage

```python
from quantum_witness import *

# 1. Build lapse fluctuation spectrum
flux_params = {'G': G, 'c': c, 'R': R, 'energy_density_psd': S_P, 'cutoff_frequency': f_c}
psd_result = build_psd_lapse_fluctuations(flux_params)

# 2. Build quantum antisymmetric spectrum
quantum_params = {'hbar': hbar, 'correlation_time': tau_c, 'noncommute_scale': alpha_nc}
antisym_result = build_antisymmetric_spectrum(psd_result, quantum_params)

# 3. Set up clock filters
clock_A = ramsey_filter(omega, T)
clock_B = hahn_echo_filter(omega, T)

# 4. Compute Q1 witness
q1_result = compute_q1_witness(clock_A, clock_B, antisym_result, psd_result['frequency'])

# 5. Verify classical limit
classical_check = verify_classical_limit({'clock_A': clock_A, 'clock_B': clock_B}, quantum_params)
```

## Verification

Run the verification script to test all components:

```bash
sage -python quantum_witness/verification.py
```

This module provides the mathematical foundation for detecting
quantum gravitational effects through precision clock networks,
enabling experimental tests of quantum spacetime structure.
