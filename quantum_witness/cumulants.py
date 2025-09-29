#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Quantum Witness: Cumulant Analysis

Implements cumulant computations K₂, K₃, K₄ from GATG Paper V.
Provides SNR analysis for higher-order gravitational signatures.

From Paper V:
K₂ = Ω² ∫ (dω/2π) |F|² S_ΔΦ
K₃ = Ω³ ∬ (dω₁dω₂/(2π)²) F₁F₂F₋₁₂ B_ΔΦ
K₄ᶜ = Ω⁴ ∭ (dω₁dω₂dω₃/(2π)³) F₁F₂F₃F₋₁₂₃ T_ΔΦ
"""

from sage.all import *

def compute_second_cumulant(filter_data, psd_data, cumulant_params=None):
    """
    Compute second cumulant K₂ (Gaussian dephasing).

    From Paper V, Equation C3:
    K₂ = Ω² ∫ (dω/2π) |F|² S_ΔΦ

    Args:
        filter_data: Dict with filter function F(ω)
        psd_data: Dict with S_ΔΦ(ω) spectrum
        cumulant_params: Optional parameters

    Returns:
        Dict with K₂ computation
    """
    if cumulant_params is None:
        cumulant_params = {}

    omega = filter_data['omega']
    F = filter_data['filter_function']
    S_Delta_Phi = psd_data['psd_lapse']
    Omega_scale = cumulant_params.get('frequency_scale', omega)

    # Power spectrum |F|²
    F_power = abs(F)**2

    # Integrand: |F|² S_ΔΦ
    integrand = F_power * S_Delta_Phi

    # K₂ integral
    K2_integral = integrate(integrand, omega, -oo, oo) / (2*pi)
    K2 = Omega_scale**2 * K2_integral

    return {
        'K2_cumulant': K2,
        'integrand': integrand,
        'frequency_scale': Omega_scale,
        'power_spectrum': F_power
    }

def compute_third_cumulant(filter_data, bispectrum_data, cumulant_params=None):
    """
    Compute third cumulant K₃ (non-Gaussian signature).

    From Paper V, Equation C3:
    K₃ = Ω³ ∬ (dω₁dω₂/(2π)²) F₁F₂F₋₁₂ B_ΔΦ

    Args:
        filter_data: Dict with filter function
        bispectrum_data: Dict with bispectrum B_ΔΦ
        cumulant_params: Optional parameters

    Returns:
        Dict with K₃ computation
    """
    if cumulant_params is None:
        cumulant_params = {}

    # Symbolic variables for double integral
    omega1 = var('omega1')
    omega2 = var('omega2')
    F = filter_data['filter_function']
    B_Delta_Phi = bispectrum_data.get('bispectrum', 0)
    Omega_scale = cumulant_params.get('frequency_scale', var('Omega'))

    # Filter products: F₁F₂F₋₁₂
    F1 = F.subs(filter_data['omega'], omega1)
    F2 = F.subs(filter_data['omega'], omega2)
    F_minus12 = F.subs(filter_data['omega'], -(omega1 + omega2))

    filter_product = F1 * F2 * F_minus12

    # K₃ double integral (symbolic)
    integrand = filter_product * B_Delta_Phi
    # Note: Actual integration requires specific bispectrum form
    K3 = Omega_scale**3 * integrand  # Symbolic form

    return {
        'K3_cumulant': K3,
        'integrand': integrand,
        'filter_product': filter_product,
        'bispectrum': B_Delta_Phi
    }

def compute_fourth_cumulant(filter_data, trispectrum_data, cumulant_params=None):
    """
    Compute fourth cumulant K₄ᶜ (connected part).

    From Paper V, Equation C3:
    K₄ᶜ = Ω⁴ ∭ (dω₁dω₂dω₃/(2π)³) F₁F₂F₃F₋₁₂₃ T_ΔΦ

    Args:
        filter_data: Dict with filter function
        trispectrum_data: Dict with trispectrum T_ΔΦ
        cumulant_params: Optional parameters

    Returns:
        Dict with K₄ᶜ computation
    """
    if cumulant_params is None:
        cumulant_params = {}

    # Symbolic variables for triple integral
    omega1, omega2, omega3 = var('omega1 omega2 omega3')
    F = filter_data['filter_function']
    T_Delta_Phi = trispectrum_data.get('trispectrum', 0)
    Omega_scale = cumulant_params.get('frequency_scale', var('Omega'))

    # Filter products: F₁F₂F₃F₋₁₂₃
    F1 = F.subs(filter_data['omega'], omega1)
    F2 = F.subs(filter_data['omega'], omega2)
    F3 = F.subs(filter_data['omega'], omega3)
    F_minus123 = F.subs(filter_data['omega'], -(omega1 + omega2 + omega3))

    filter_product = F1 * F2 * F3 * F_minus123

    # K₄ᶜ triple integral (symbolic)
    integrand = filter_product * T_Delta_Phi
    K4c = Omega_scale**4 * integrand  # Symbolic form

    return {
        'K4c_cumulant': K4c,
        'integrand': integrand,
        'filter_product': filter_product,
        'trispectrum': T_Delta_Phi
    }

def compute_snr_lines(cumulant_results, measurement_params):
    """
    Compute SNR lines for cumulant detection.

    From Paper V, Equation X1:
    Var K̂₂ = 2σ_φ⁴/M, Var K̂₃ = 6σ_φ⁶/M, Var K̂₄ = 24σ_φ⁸/M
    SNR₂ = √M K₂/(√2 σ_φ²), SNR₃ = √M K₃/(√6 σ_φ³), SNR₄ = √M K₄/(√24 σ_φ⁴)

    Args:
        cumulant_results: Dict with computed cumulants
        measurement_params: Dict with M, σ_φ

    Returns:
        Dict with SNR computations
    """
    M = measurement_params['num_measurements']
    sigma_phi = measurement_params['per_shot_phase_noise']

    K2 = cumulant_results.get('K2_cumulant', 0)
    K3 = cumulant_results.get('K3_cumulant', 0)
    K4c = cumulant_results.get('K4c_cumulant', 0)

    # Variances from shot noise
    Var_K2 = 2 * sigma_phi**4 / M
    Var_K3 = 6 * sigma_phi**6 / M
    Var_K4 = 24 * sigma_phi**8 / M

    # SNR calculations
    SNR2 = sqrt(M) * K2 / (sqrt(2) * sigma_phi**2)
    SNR3 = sqrt(M) * K3 / (sqrt(6) * sigma_phi**3)
    SNR4 = sqrt(M) * K4c / (sqrt(24) * sigma_phi**4)

    return {
        'SNR2': SNR2,
        'SNR3': SNR3,
        'SNR4': SNR4,
        'variances': {
            'Var_K2': Var_K2,
            'Var_K3': Var_K3,
            'Var_K4': Var_K4
        },
        'measurement_count': M,
        'phase_noise': sigma_phi
    }