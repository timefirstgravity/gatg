#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Dephasing Observatory: Allan Variance Bridge

Allan variance computations and PSD conversions.
Provides bridge between frequency-domain gravitational signatures
and time-domain clock stability measurements.

From GATG Paper V: Allan variance bridge for clock network analysis.
"""

from sage.all import *

def compute_allan_deviation(psd_data, averaging_times, allan_params=None):
    """
    Compute Allan deviation from power spectral density.

    From Paper V, Equation N1:
    σ_y²(τ_m) = 2 ∫₀^∞ df S^(1)_y(f) sin⁴(π f τ_m) / (π f τ_m)²

    Args:
        psd_data: Dict with one-sided PSD S^(1)_y(f)
        averaging_times: List of averaging times τ_m
        allan_params: Optional computation parameters

    Returns:
        Dict with Allan deviation values and analysis
    """
    if allan_params is None:
        allan_params = {'method': 'symbolic'}

    # Extract PSD
    f = psd_data['frequency']
    S_y_one_sided = psd_data['psd_one_sided']

    # Allan deviation computation
    allan_variances = {}
    allan_deviations = {}

    for tau_m in averaging_times:
        # Allan variance integrand: 2 S^(1)_y(f) sin⁴(π f τ_m) / (π f τ_m)²
        sin_factor = sin(pi * f * tau_m)**4
        normalization = (pi * f * tau_m)**2
        integrand = 2 * S_y_one_sided * sin_factor / normalization

        # Integration from 0 to infinity
        if allan_params['method'] == 'symbolic':
            allan_variance = integrate(integrand, f, 0, oo)
        else:
            raise ValueError(f"Unknown Allan method: {allan_params['method']}")

        allan_variances[tau_m] = allan_variance
        allan_deviations[tau_m] = sqrt(allan_variance)

    return {
        'allan_variances': allan_variances,
        'allan_deviations': allan_deviations,
        'averaging_times': averaging_times,
        'psd_data': psd_data
    }

def convert_psd_to_allan(two_sided_psd, frequency_var):
    """
    Convert two-sided PSD to Allan variance format.

    Converts S(f) (two-sided) to S^(1)(f) (one-sided) for Allan computation.

    Args:
        two_sided_psd: Two-sided power spectral density
        frequency_var: Frequency variable

    Returns:
        Dict with one-sided PSD for Allan variance computation
    """
    f = frequency_var

    # One-sided PSD: S^(1)(f) = 2 S(f) for f > 0
    psd_one_sided = 2 * two_sided_psd

    return {
        'psd_one_sided': psd_one_sided,
        'frequency': f,
        'conversion': 'two_sided_to_one_sided',
        'factor': 2
    }

def verify_allan_scaling(allan_result, expected_scaling):
    """
    Verify Allan deviation scaling with averaging time.

    Different noise types have characteristic Allan variance scaling:
    - White noise: σ_y(τ) ∝ 1/√τ
    - Flicker noise: σ_y(τ) ∝ constant
    - Random walk: σ_y(τ) ∝ √τ

    Args:
        allan_result: Result from compute_allan_deviation
        expected_scaling: Dict with expected scaling parameters

    Returns:
        Dict with scaling verification
    """
    allan_devs = allan_result['allan_deviations']
    tau_list = allan_result['averaging_times']

    # Extract scaling exponent from Allan deviation vs τ
    # σ_y(τ) = A τ^α ⟹ log σ_y = log A + α log τ

    scaling_exponents = {}
    if len(tau_list) >= 2:
        # Compute scaling between consecutive points
        for i in range(len(tau_list) - 1):
            tau1, tau2 = tau_list[i], tau_list[i + 1]
            sigma1, sigma2 = allan_devs[tau1], allan_devs[tau2]

            # Scaling exponent: α = log(σ₂/σ₁) / log(τ₂/τ₁)
            scaling_exp = log(sigma2/sigma1) / log(tau2/tau1)
            scaling_exponents[f'tau_{tau1}_to_{tau2}'] = scaling_exp

    # Compare with expected scaling
    expected_exp = expected_scaling.get('exponent', 0)
    noise_type = expected_scaling.get('type', 'unknown')

    scaling_verified = {}
    for key, actual_exp in scaling_exponents.items():
        # Check if actual scaling matches expected (within tolerance)
        tolerance = expected_scaling.get('tolerance', 0.1)
        scaling_verified[key] = abs(actual_exp - expected_exp) < tolerance

    return {
        'scaling_exponents': scaling_exponents,
        'expected_exponent': expected_exp,
        'noise_type': noise_type,
        'scaling_verified': scaling_verified,
        'overall_scaling_correct': all(scaling_verified.values()) if scaling_verified else None
    }

def analyze_noise_type(allan_result, tau_range):
    """
    Analyze noise type from Allan deviation scaling.

    Identifies dominant noise processes from Allan variance behavior.

    Args:
        allan_result: Result from compute_allan_deviation
        tau_range: Range of averaging times for analysis

    Returns:
        Dict with noise type identification
    """
    allan_devs = allan_result['allan_deviations']

    # Standard noise identification by Allan variance scaling
    noise_classifications = {
        'white_phase': {'exponent': -1, 'name': 'White Phase Noise'},
        'flicker_phase': {'exponent': -0.5, 'name': 'Flicker Phase Noise'},
        'white_frequency': {'exponent': 0, 'name': 'White Frequency Noise'},
        'flicker_frequency': {'exponent': 0.5, 'name': 'Flicker Frequency Noise'},
        'random_walk': {'exponent': 1, 'name': 'Random Walk Frequency Noise'}
    }

    # Compute average scaling exponent over specified range
    tau_in_range = [tau for tau in allan_result['averaging_times'] if tau_range[0] <= tau <= tau_range[1]]

    if len(tau_in_range) >= 2:
        # Fit Allan deviation to power law: σ_y(τ) = A τ^α
        # Use first and last points in range
        tau_start, tau_end = tau_in_range[0], tau_in_range[-1]
        sigma_start, sigma_end = allan_devs[tau_start], allan_devs[tau_end]

        # Scaling exponent
        measured_exponent = log(sigma_end/sigma_start) / log(tau_end/tau_start)

        # Find closest noise type
        closest_noise = None
        min_difference = float('inf')
        for noise_key, noise_info in noise_classifications.items():
            diff = abs(measured_exponent - noise_info['exponent'])
            if diff < min_difference:
                min_difference = diff
                closest_noise = noise_key

        noise_identification = {
            'measured_exponent': measured_exponent,
            'closest_noise_type': closest_noise,
            'noise_name': noise_classifications[closest_noise]['name'],
            'exponent_difference': min_difference
        }
    else:
        noise_identification = {
            'measured_exponent': None,
            'closest_noise_type': 'insufficient_data',
            'noise_name': 'Cannot determine from data',
            'exponent_difference': None
        }

    return {
        'noise_identification': noise_identification,
        'noise_classifications': noise_classifications,
        'analysis_range': tau_range,
        'confidence': 'high' if min_difference < 0.1 else 'low'
    }

def build_allan_bridge_to_psd(allan_measurements, bridge_params=None):
    """
    Build bridge from Allan variance measurements back to PSD estimation.

    Inverse operation: estimate PSD from Allan variance data.

    Args:
        allan_measurements: Dict with measured Allan variance data
        bridge_params: Optional bridge computation parameters

    Returns:
        Dict with estimated PSD and bridge analysis
    """
    if bridge_params is None:
        bridge_params = {'model': 'power_law'}

    tau_values = allan_measurements['averaging_times']
    sigma_values = allan_measurements['allan_deviations']

    # Model Allan variance as power law: σ²(τ) = A τ^β
    tau_var = var('tau', domain='positive')
    f_var = var('f', domain='positive')

    # For power law Allan variance, corresponding PSD has specific form
    if bridge_params['model'] == 'power_law':
        # Fit power law to Allan data (symbolic)
        # σ²(τ) = A τ^β ⟺ S_y(f) = B f^γ where γ = -(β + 1)

        # Use first and last data points to estimate power law
        if len(tau_values) >= 2:
            tau1, tau2 = tau_values[0], tau_values[-1]
            sigma1, sigma2 = sigma_values[tau1], sigma_values[tau2]

            # Power law exponent
            beta = log(sigma2**2 / sigma1**2) / log(tau2 / tau1)
            A = sigma1**2 / tau1**beta

            # Corresponding PSD exponent
            gamma = -(beta + 1)

            # Estimated PSD: S_y(f) = B f^γ
            # Bridge coefficient B determined by Allan-PSD relationship
            B = A * (2 * log(2)) / (pi**2) if beta == -1 else A  # Symbolic form

            estimated_psd = B * f_var**gamma

        else:
            estimated_psd = var('estimated_psd')
            beta = var('beta')
            gamma = var('gamma')

    else:
        raise ValueError(f"Unknown bridge model: {bridge_params['model']}")

    return {
        'estimated_psd': estimated_psd,
        'power_law_exponent': gamma,
        'allan_exponent': beta,
        'bridge_model': bridge_params['model'],
        'frequency_variable': f_var
    }