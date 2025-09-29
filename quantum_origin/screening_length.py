#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Quantum Origin: Screening Length Extraction

Extract the gravitational screening length ξ from CTP kernel spectrum.
The screening length determines the range of gravitational correlations
and must satisfy ξ ≳ 10¹¹ m for consistency with solar system tests.

From Quantum Origin paper:
- Screened Poisson: (-∇² + m²)Φ = (4πG/c⁴)ρ
- Screening length: ξ = 1/m
- Green's function tail: G(r) ~ exp(-r/ξ)/r
"""

from sage.all import *
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def extract_screening_length(kernel_data, extraction_params=None):
    """
    Extract screening length ξ from kernel properties.

    Methods:
    1. Direct from screening mass: ξ = 1/m
    2. From Green's function decay: G(r) ~ exp(-r/ξ)
    3. From spectral gap: infrared regularization scale

    Args:
        kernel_data: CTP kernel with spectral properties
        extraction_params: Method and parameters

    Returns:
        Dict with:
            'screening_length': ξ in meters
            'screening_mass': m = 1/ξ
            'extraction_method': Method used
            'observational_bound_satisfied': ξ ≥ 10¹¹ m
    """
    if extraction_params is None:
        extraction_params = {'method': 'auto'}

    kernel_type = kernel_data['kernel_type']
    method = extraction_params['method']

    xi_min = 1e11  # Minimum 10^11 meters from solar system tests

    if kernel_type == 'screened':
        # Direct extraction from screened kernel
        xi = kernel_data['screening_length']
        m = kernel_data['screening_mass']
        method_used = 'direct_from_screening'

    elif kernel_type == 'quasilocal':
        # Extract from correlation length
        ell_c = kernel_data['correlation_length']
        # Screening length related to correlation length
        # Factor depends on specific kernel form
        xi = ell_c  # Simplified - actual factor depends on kernel
        m = 1/xi if xi > 0 else float('inf')
        method_used = 'from_correlation_length'

    elif kernel_type == 'local':
        # No screening for local kernel
        xi = float('inf')  # Infinite screening length
        m = 0  # Zero screening mass
        method_used = 'no_screening_local_kernel'

    else:
        # Try to extract from spectrum
        if 'green_spectrum' in kernel_data:
            # From momentum space: G̃(k) = A/(k² + m²)
            # The pole at k² = -m² gives screening mass
            spectrum = kernel_data['green_spectrum']
            k = var('k')

            # Find where spectrum has pole (simplified)
            # In practice, would solve for pole location
            m = kernel_data.get('screening_mass', 0)
            xi = 1/m if m > 0 else float('inf')
            method_used = 'from_spectrum_pole'
        else:
            xi = None
            m = None
            method_used = 'extraction_failed'

    # Check observational bound
    if xi is not None:
        bound_satisfied = xi >= xi_min if xi != float('inf') else True
    else:
        bound_satisfied = None

    return {
        'screening_length': xi,
        'screening_mass': m,
        'extraction_method': method_used,
        'observational_bound_satisfied': bound_satisfied,
        'minimum_required': xi_min,
        'units': 'meters',
        'physical_interpretation': 'gravitational_correlation_range'
    }

def compute_correlation_function(kernel_data, spatial_points, correlation_params=None):
    """
    Compute spatial correlation function from kernel.

    The correlation function C(r) = ⟨Φ(x)Φ(x+r)⟩ determines
    how gravitational effects correlate over distance.

    Args:
        kernel_data: CTP kernel data
        spatial_points: List of distances r to evaluate
        correlation_params: Optional parameters

    Returns:
        Dict with correlation function values
    """
    if correlation_params is None:
        correlation_params = {}

    kernel_type = kernel_data['kernel_type']

    correlations = {}

    if kernel_type == 'local':
        # Delta function correlation
        for r in spatial_points:
            correlations[r] = 1.0 if r == 0 else 0.0

    elif kernel_type == 'quasilocal':
        # Exponential or Gaussian decay
        kernel_func = kernel_data['kernel_function']
        r_var = var('r')

        for r_val in spatial_points:
            # Evaluate kernel at distance r
            corr_value = kernel_func.subs(r_var == r_val)
            correlations[r_val] = corr_value

    elif kernel_type == 'screened':
        # Yukawa-type correlation
        green_func = kernel_data['green_function']
        r_var = var('r')

        for r_val in spatial_points:
            # Green's function gives correlation
            corr_value = green_func.subs(r_var == r_val)
            correlations[r_val] = corr_value

    # Analyze decay
    if len(spatial_points) > 1:
        # Fit exponential decay to extract length scale
        r_vals = sorted(spatial_points)
        if kernel_type != 'local':
            # Check for exponential decay pattern
            first_nonzero = correlations[r_vals[0]] if r_vals[0] > 0 else correlations[r_vals[1]]
            last = correlations[r_vals[-1]]

            if first_nonzero != 0 and last != 0:
                # Decay length from ratio
                decay_length = (r_vals[-1] - r_vals[0]) / log(abs(first_nonzero/last))
            else:
                decay_length = None
        else:
            decay_length = 0  # Local kernel
    else:
        decay_length = None

    return {
        'correlation_function': correlations,
        'spatial_points': spatial_points,
        'decay_length': decay_length,
        'kernel_type': kernel_type,
        'long_range_behavior': 'exponential' if kernel_type == 'screened' else 'power_law'
    }

def analyze_greens_function_tail(kernel_data, large_r_limit=1e6, tail_params=None):
    """
    Analyze Green's function tail behavior at large distances.

    The tail G(r→∞) determines:
    - Screening length from exponential decay
    - Power law exponents for unscreened cases
    - Observational constraints on modifications

    Args:
        kernel_data: Kernel with Green's function
        large_r_limit: Distance for asymptotic analysis
        tail_params: Optional analysis parameters

    Returns:
        Dict with tail behavior analysis
    """
    if tail_params is None:
        tail_params = {}

    kernel_type = kernel_data['kernel_type']

    if kernel_type == 'screened':
        # Yukawa tail: G(r) ~ exp(-mr)/r
        m = kernel_data['screening_mass']
        xi = kernel_data['screening_length']

        r = var('r', domain='positive')
        green_func = kernel_data['green_function']

        # Leading asymptotic behavior
        # G(r→∞) ~ A exp(-r/ξ)/r
        amplitude_prefactor = kernel_data['spectral_properties'].get('max_eigenvalue', 1)

        tail_form = amplitude_prefactor * exp(-r/xi) / r

        # Decay rate
        exponential_decay_rate = 1/xi

        # Distance where G falls to 1/e of Newtonian
        decay_distance = xi

        tail_analysis = {
            'asymptotic_form': str(tail_form),
            'decay_type': 'exponential',
            'decay_rate': exponential_decay_rate,
            'decay_distance': decay_distance,
            'screening_length': xi
        }

    elif kernel_type == 'quasilocal':
        # Depends on specific kernel form
        ell_c = kernel_data['correlation_length']

        if kernel_data.get('decay_type') == 'exponential':
            tail_analysis = {
                'asymptotic_form': 'exp(-r/ell_c)/r',
                'decay_type': 'exponential',
                'decay_rate': 1/ell_c,
                'decay_distance': ell_c,
                'correlation_length': ell_c
            }
        elif kernel_data.get('decay_type') == 'gaussian':
            tail_analysis = {
                'asymptotic_form': 'exp(-r²/2ell_c²)',
                'decay_type': 'gaussian',
                'decay_rate': 1/ell_c**2,
                'decay_distance': ell_c,
                'correlation_length': ell_c
            }
        else:
            tail_analysis = {'decay_type': 'unknown'}

    elif kernel_type == 'local':
        # No tail for local kernel
        tail_analysis = {
            'asymptotic_form': 'delta(r)',
            'decay_type': 'none',
            'decay_rate': float('inf'),
            'decay_distance': 0
        }

    else:
        tail_analysis = {'decay_type': 'undetermined'}

    # Evaluate at large r
    if 'green_function' in kernel_data and kernel_type != 'local':
        r_var = var('r')
        green_at_large_r = kernel_data['green_function'].subs(r_var == large_r_limit)
        tail_analysis['value_at_large_r'] = green_at_large_r
        tail_analysis['evaluation_distance'] = large_r_limit

    return tail_analysis

def verify_screening_bounds(screening_data, observational_constraints=None):
    """
    Verify screening length against observational bounds.

    Constraints:
    1. Solar system: ξ > 10¹¹ m (Cassini bound)
    2. Galactic: ξ > 10¹⁶ m (rotation curves)
    3. Cosmological: ξ < Hubble radius

    Args:
        screening_data: Result from extract_screening_length
        observational_constraints: Dict with specific bounds

    Returns:
        Dict with constraint verification
    """
    if observational_constraints is None:
        observational_constraints = {
            'solar_system_min': 1e11,      # 10^11 m (Cassini)
            'galactic_min': 1e16,           # 10^16 m (optional stronger)
            'cosmological_max': 1e26        # Hubble radius
        }

    xi = screening_data['screening_length']

    if xi is None:
        return {
            'all_bounds_satisfied': None,
            'message': 'No screening length to verify'
        }

    # Check each bound
    solar_satisfied = xi >= observational_constraints['solar_system_min']

    if xi == float('inf'):
        # Infinite screening length (no screening)
        galactic_satisfied = True
        cosmological_satisfied = True
        interpretation = 'No screening (standard GR)'
    else:
        galactic_satisfied = xi >= observational_constraints['galactic_min']
        cosmological_satisfied = xi <= observational_constraints['cosmological_max']

        if xi < 1e11:
            interpretation = 'Too short - ruled out by solar system tests'
        elif xi < 1e16:
            interpretation = 'Solar system compatible, may affect galactic dynamics'
        elif xi < 1e26:
            interpretation = 'Compatible with all current tests'
        else:
            interpretation = 'Very long range - essentially unscreened'

    all_satisfied = solar_satisfied and cosmological_satisfied

    return {
        'all_bounds_satisfied': all_satisfied,
        'solar_system_bound_satisfied': solar_satisfied,
        'galactic_bound_satisfied': galactic_satisfied,
        'cosmological_bound_satisfied': cosmological_satisfied,
        'screening_length': xi,
        'bounds_used': observational_constraints,
        'physical_interpretation': interpretation,
        'units': 'meters'
    }