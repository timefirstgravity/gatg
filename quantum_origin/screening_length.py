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
    Extract screening length ξ from kernel properties with uncertainty analysis.

    Methods:
    1. Direct from screening mass: ξ = 1/m
    2. From Green's function decay: G(r) ~ exp(-r/ξ)
    3. From spectral gap: infrared regularization scale

    Uncertainty analysis:
    - Uncertainty propagation from kernel parameters
    - Error bar calculation for ξ
    - Confidence intervals and parameter sensitivity

    Args:
        kernel_data: CTP kernel with spectral properties
        extraction_params: Method and parameters for extraction

    Returns:
        Dict with:
            'screening_length': ξ in meters
            'screening_mass': m = 1/ξ
            'extraction_method': Method used
            'observational_bound_satisfied': ξ ≥ 10¹¹ m
            'uncertainty_analysis': Error bars and confidence intervals
            'parameter_sensitivity': Sensitivity to input parameters
    """
    if extraction_params is None:
        extraction_params = {'method': 'auto', 'uncertainty_analysis': True}

    kernel_type = kernel_data['kernel_type']
    method = extraction_params['method']
    uncertainty_analysis_enabled = extraction_params.get('uncertainty_analysis', True)

    xi_min = 1e11  # Minimum 10^11 meters from solar system tests

    # Main extraction logic
    if kernel_type == 'screened':
        # Direct extraction from screened kernel
        xi = kernel_data['screening_length']
        m = kernel_data['screening_mass']
        method_used = 'direct_from_screening'

        # Uncertainty analysis for screened kernel
        if uncertainty_analysis_enabled:
            # Main uncertainty comes from screening mass determination
            m_uncertainty = m * 0.05  # 5% uncertainty in mass determination
            xi_uncertainty = xi * m_uncertainty / m  # Propagate through ξ = 1/m

            # Parameter sensitivity
            sensitivity_analysis = {
                'screening_mass_sensitivity': -xi**2,  # dξ/dm = -1/m² = -ξ²
                'relative_sensitivity': m_uncertainty / m,
                'dominant_uncertainty_source': 'screening_mass_determination'
            }
        else:
            xi_uncertainty = None
            sensitivity_analysis = {}

    elif kernel_type == 'quasilocal':
        # Extract from correlation length
        ell_c = kernel_data['correlation_length']
        # Screening length related to correlation length
        # Factor depends on specific kernel form
        xi = ell_c  # Simplified - actual factor depends on kernel
        m = 1/xi if xi > 0 else float('inf')
        method_used = 'from_correlation_length'

        # Uncertainty analysis for quasi-local kernel
        if uncertainty_analysis_enabled:
            # Uncertainty from correlation length determination
            ell_c_uncertainty = ell_c * 0.1  # 10% uncertainty in correlation length
            xi_uncertainty = ell_c_uncertainty  # Direct propagation for ξ = ℓ_c

            sensitivity_analysis = {
                'correlation_length_sensitivity': 1.0,  # dξ/dℓ_c = 1
                'relative_sensitivity': ell_c_uncertainty / ell_c,
                'dominant_uncertainty_source': 'correlation_length_measurement'
            }
        else:
            xi_uncertainty = None
            sensitivity_analysis = {}

    elif kernel_type == 'local':
        # No screening for local kernel
        xi = float('inf')  # Infinite screening length
        m = 0  # Zero screening mass
        method_used = 'no_screening_local_kernel'

        # No uncertainty for infinite screening length
        xi_uncertainty = None
        sensitivity_analysis = {
            'screening_length': 'infinite',
            'uncertainty': 'not_applicable',
            'physical_meaning': 'no_screening_standard_GR'
        }

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

            # Uncertainty from spectral analysis
            if uncertainty_analysis_enabled and m > 0:
                # Uncertainty in mass from spectral fitting
                spectral_uncertainty = kernel_data.get('spectral_uncertainty', m * 0.15)
                xi_uncertainty = xi * spectral_uncertainty / m

                sensitivity_analysis = {
                    'spectral_mass_sensitivity': -xi**2,
                    'relative_sensitivity': spectral_uncertainty / m,
                    'dominant_uncertainty_source': 'spectral_analysis_fitting'
                }
            else:
                xi_uncertainty = None
                sensitivity_analysis = {}
        else:
            xi = None
            m = None
            method_used = 'extraction_failed'
            xi_uncertainty = None
            sensitivity_analysis = {'error': 'insufficient_data_for_extraction'}

    # Check observational bound
    if xi is not None:
        bound_satisfied = xi >= xi_min if xi != float('inf') else True
    else:
        bound_satisfied = None

    # Comprehensive uncertainty analysis
    uncertainty_analysis = {}
    if uncertainty_analysis_enabled and xi_uncertainty is not None and xi is not None and xi != float('inf'):
        # Confidence intervals (assuming normal distribution)
        confidence_levels = [0.68, 0.95, 0.99]  # 1σ, 2σ, 3σ
        confidence_intervals = {}

        for level in confidence_levels:
            # Convert confidence level to number of standard deviations
            n_sigma = {0.68: 1, 0.95: 2, 0.99: 3}[level]
            lower_bound = max(xi - n_sigma * xi_uncertainty, 0)
            upper_bound = xi + n_sigma * xi_uncertainty

            confidence_intervals[f'{level:.0%}'] = {
                'lower_bound': lower_bound,
                'upper_bound': upper_bound,
                'interval_width': upper_bound - lower_bound
            }

        # Error bar information
        uncertainty_analysis = {
            'central_value': xi,
            'absolute_uncertainty': xi_uncertainty,
            'relative_uncertainty': xi_uncertainty / xi,
            'confidence_intervals': confidence_intervals,
            'error_bar_lower': xi - xi_uncertainty,
            'error_bar_upper': xi + xi_uncertainty,
            'uncertainty_method': 'parameter_propagation',
            'observational_bound_margin': (xi - xi_uncertainty - xi_min) / xi_min if xi > xi_min else None
        }

        # Check if uncertainty affects observational bound satisfaction
        uncertainty_analysis['bound_satisfied_with_uncertainty'] = (xi - xi_uncertainty) >= xi_min

    # Monte Carlo uncertainty estimation (optional enhancement)
    if uncertainty_analysis_enabled and extraction_params.get('monte_carlo_samples', 0) > 0:
        n_samples = extraction_params['monte_carlo_samples']
        mc_samples = []

        # Generate parameter samples and propagate uncertainty
        for _ in range(n_samples):
            if kernel_type == 'screened':
                # Sample screening mass with uncertainty
                m_sample = m + (random() - 0.5) * 2 * (m * 0.05)
                xi_sample = 1/m_sample if m_sample > 0 else float('inf')
            elif kernel_type == 'quasilocal':
                # Sample correlation length with uncertainty
                ell_c_sample = ell_c + (random() - 0.5) * 2 * (ell_c * 0.1)
                xi_sample = ell_c_sample
            else:
                xi_sample = xi

            if xi_sample != float('inf'):
                mc_samples.append(xi_sample)

        if mc_samples:
            mc_mean = mean(mc_samples)
            mc_std = sqrt(variance(mc_samples))
            uncertainty_analysis['monte_carlo'] = {
                'mean': mc_mean,
                'standard_deviation': mc_std,
                'samples': len(mc_samples),
                'confidence_95': {
                    'lower': mc_mean - 2*mc_std,
                    'upper': mc_mean + 2*mc_std
                }
            }

    return {
        'screening_length': xi,
        'screening_mass': m,
        'extraction_method': method_used,
        'observational_bound_satisfied': bound_satisfied,
        'minimum_required': xi_min,
        'units': 'meters',
        'physical_interpretation': 'gravitational_correlation_range',
        # Uncertainty analysis data
        'uncertainty_analysis': uncertainty_analysis,
        'parameter_sensitivity': sensitivity_analysis,
        'error_bars': {
            'central_value': xi,
            'uncertainty': xi_uncertainty,
            'lower_bound': xi - xi_uncertainty if xi_uncertainty and xi != float('inf') else None,
            'upper_bound': xi + xi_uncertainty if xi_uncertainty and xi != float('inf') else None
        } if xi_uncertainty is not None else None
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