# Dephasing Observatory Module

Implements gravitational dephasing predictions for clock networks based on GATG Quantum Origin theory. This module provides rigorous mathematical computations of capacity Ξ, visibility V, and linewidth Δν for gravitational decoherence effects in precision clock experiments.

## Core Physics

The module quantifies how gravitational field fluctuations cause decoherence in separated quantum clocks through three key observables:

1. **Capacity** Ξ: Accumulated phase noise from lapse fluctuations
2. **Visibility** V: Coherence loss in clock comparisons
3. **Linewidth** Δν: Frequency broadening from gravitational noise

## Mathematical Foundation

### Key Equations

**Capacity**: `Ξ = ∫ (dω/2π) |W_T(ω)|² S_Φ(ω)`

**Visibility**: `V = exp[-½ω²Ξ]`

**Linewidth**: `Δν = (ω/2πT)√Ξ`

Where:
- `S_Φ(ω)` is the lapse fluctuation power spectral density
- `W_T(ω)` is the observation window function
- `ω` is the clock angular frequency
- `T` is the observation time

## Module Components

### `capacity.py`
Core capacity computations:
- `compute_capacity_xi()` - Capacity Ξ from lapse PSD and window function
- `compute_windowed_capacity()` - Capacity with rectangular/Hann/Gaussian windows
- `verify_capacity_properties()` - Mathematical and physical verification
- `compute_correlation_length_xi()` - Extract screening correlation length

### `visibility.py`
Visibility and coherence analysis:
- `compute_visibility()` - Visibility V from capacity and clock frequency
- `compute_linewidth()` - Linewidth Δν from gravitational dephasing
- `compute_dephasing_rate()` - Dephasing rate Γ characterization
- `verify_exponential_scaling()` - Verify V = exp[-½ω²Ξ] properties
- `compute_visibility_contours()` - Parameter space mapping

### `snr_analysis.py`
Signal-to-noise and parameter estimation:
- `compute_fisher_matrix()` - Fisher information for parameter estimation
- `compute_parameter_errors()` - Error estimates from Fisher matrix
- `build_snr_isolines()` - SNR contours in parameter space
- `verify_cramer_rao_bound()` - Fundamental limit verification
- `optimize_measurement_strategy()` - Experimental design optimization

### `allan_variance.py`
Clock stability metrics bridge:
- `compute_allan_deviation()` - Allan deviation from PSD
- `convert_psd_to_allan()` - Two-sided to one-sided PSD conversion
- `verify_allan_scaling()` - Noise type verification
- `analyze_noise_type()` - Identify dominant noise processes
- `build_allan_bridge_to_psd()` - Inverse: estimate PSD from Allan data

### `scenarios.py`
Realistic experimental scenarios:
- `leo_ground_scenario()` - LEO satellite to ground clock comparison
- `ground_ground_scenario()` - Separated ground-based clocks
- `intercontinental_scenario()` - Global clock network analysis
- `verify_physical_constraints()` - Causality and detectability checks

## Physical Scenarios

### LEO-Ground (400 km altitude)
- Baseline: ~400 km
- Earth's gravitational field variations
- Detectable with optical clocks
- Observation time: ~1000 s

### Ground-Ground (1-1000 km)
- Geological density fluctuations dominate
- Distance scaling: Ξ ∝ 1/L²
- Current clock stability sufficient
- Observation time: ~1 hour

### Intercontinental (10,000+ km)
- Global field effects (tidal, atmospheric, tectonic)
- Earth curvature corrections included
- Long observation times (24+ hours)
- Tests GR at planetary scales

## Key Features

✓ **Rigorous symbolic computation** with SageMath
✓ **No synthetic data** - all calculations from first principles
✓ **Physical regularization** - Lorentzian cutoffs prevent divergences
✓ **Multiple window functions** - rectangular, Hann, Gaussian
✓ **Complete verification suite** - dimensional, scaling, causality checks
✓ **Realistic parameters** - uses actual G, c, Earth mass values

## Usage Example

```python
from dephasing_observatory import *

# 1. Build lapse PSD from flux law
flux_params = {
    'G': 6.67430e-11,  # SI units
    'c': 299792458,
    'R': 1000,  # 1 km baseline
    'energy_density_psd': S_P,
    'cutoff_frequency': 1e-3  # 1 mHz
}
psd_data = build_psd_lapse_fluctuations(flux_params)

# 2. Compute windowed capacity
window_params = {'T': 1000}  # 1000 s observation
capacity_result = compute_windowed_capacity(psd_data, 'rectangular', window_params)

# 3. Calculate visibility loss
omega = 2*pi * 1e15  # Optical clock frequency
visibility_result = compute_visibility(capacity_result['capacity_xi'], omega)

# 4. Determine linewidth broadening
linewidth_result = compute_linewidth(capacity_result['capacity_xi'], omega, T_obs=1000)

# 5. Run physical scenario
satellite_params = {'altitude': 400000, 'clock_frequency': omega}
ground_params = {'clock_frequency': omega}
leo_result = leo_ground_scenario(satellite_params, ground_params)
```

## Experimental Predictions

The module provides testable predictions for:
- **Visibility loss** in clock comparisons
- **Frequency broadening** from gravitational noise
- **Distance scaling** of decoherence effects
- **Optimal measurement strategies** for detection

## Verification

Run the comprehensive verification script:
```bash
sage -python dephasing_observatory/verification.py
```

This module enables precision tests of gravitational decoherence in quantum clock networks, providing the mathematical foundation for experimental validation of GATG predictions about the quantum nature of spacetime.