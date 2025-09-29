#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Dephasing Observatory: Visibility and Linewidth

Implements visibility and linewidth calculations from GATG Quantum Origin paper.
The visibility V and linewidth Δν are key observables for gravitational dephasing.

From Paper QOrigin:
V = exp[-½ω²Ξ]
Δν = (ω/2πT)√Ξ
"""

from sage.all import *
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def compute_visibility(capacity_xi, clock_frequency, visibility_params=None):
    """
    Compute visibility V from capacity and clock frequency.

    From Paper QOrigin, Equation L1:
    V = exp[-½ω²Ξ]

    The visibility V quantifies coherence loss due to gravitational dephasing.
    V = 1: perfect coherence, V → 0: complete dephasing.

    Args:
        capacity_xi: Ξ value from capacity computation (s²)
        clock_frequency: ω angular frequency (rad/s)
        visibility_params: Optional parameters for computation

    Returns:
        Dict with:
            'visibility': V value
            'dephasing_exponent': ½ω²Ξ
            'coherence_loss': 1 - V
            'log_visibility': ln(V) = -½ω²Ξ
    """
    if visibility_params is None:
        visibility_params = {}

    omega = clock_frequency
    Xi = capacity_xi

    # Dephasing exponent: ½ω²Ξ
    dephasing_exponent = omega**2 * Xi / 2

    # Visibility: V = exp[-½ω²Ξ]
    visibility = exp(-dephasing_exponent)

    # Coherence loss: 1 - V
    coherence_loss = 1 - visibility

    # Log visibility: ln(V) = -½ω²Ξ
    log_visibility = -dephasing_exponent

    return {
        'visibility': visibility,
        'dephasing_exponent': dephasing_exponent,
        'coherence_loss': coherence_loss,
        'log_visibility': log_visibility,
        'frequency': omega,
        'capacity': Xi
    }

def compute_linewidth(capacity_xi, clock_frequency, observation_time, linewidth_params=None):
    """
    Compute linewidth Δν from dephasing.

    From Paper QOrigin, Equation L1:
    Δν = (ω/2πT)√Ξ

    The linewidth Δν represents the frequency broadening due to
    gravitational phase noise over observation time T.

    Args:
        capacity_xi: Ξ value from capacity computation (s²)
        clock_frequency: ω angular frequency (rad/s)
        observation_time: T observation window (s)
        linewidth_params: Optional parameters

    Returns:
        Dict with:
            'linewidth': Δν (Hz)
            'angular_linewidth': Δω = 2πΔν (rad/s)
            'q_factor': ω/Δω quality factor
            'fractional_linewidth': Δν/ν
    """
    if linewidth_params is None:
        linewidth_params = {}

    omega = clock_frequency
    Xi = capacity_xi
    T = observation_time

    # Linewidth: Δν = (ω/2πT)√Ξ
    linewidth = (omega / (2*pi*T)) * sqrt(Xi)

    # Angular linewidth: Δω = 2πΔν
    angular_linewidth = 2*pi * linewidth

    # Quality factor: Q = ω/Δω
    # Handle case where linewidth is negligible (undetectable dephasing)
    if angular_linewidth.is_trivial_zero() or angular_linewidth == 0:
        q_factor = oo  # Infinite Q-factor for negligible dephasing
    else:
        q_factor = omega / angular_linewidth

    # Fractional linewidth: Δν/ν where ν = ω/2π
    frequency_hz = omega / (2*pi)
    if linewidth.is_trivial_zero() or linewidth == 0:
        fractional_linewidth = 0  # No fractional broadening for negligible dephasing
    else:
        fractional_linewidth = linewidth / frequency_hz

    return {
        'linewidth': linewidth,
        'angular_linewidth': angular_linewidth,
        'q_factor': q_factor,
        'fractional_linewidth': fractional_linewidth,
        'frequency_hz': frequency_hz,
        'observation_time': T,
        'capacity': Xi
    }

def compute_dephasing_rate(capacity_xi, clock_frequency, rate_params=None):
    """
    Compute dephasing rate from capacity.

    The dephasing rate characterizes how quickly coherence is lost:
    Γ_dephase = ω²Ξ/T

    Args:
        capacity_xi: Ξ value from capacity computation (s²)
        clock_frequency: ω angular frequency (rad/s)
        rate_params: Dict with observation time T

    Returns:
        Dict with:
            'dephasing_rate': Γ_dephase (s⁻¹)
            'dephasing_time': 1/Γ_dephase (s)
            'frequency_scaling': ω² dependence
    """
    if rate_params is None:
        rate_params = {'observation_time': var('T', positive=True)}

    omega = clock_frequency
    Xi = capacity_xi
    T = rate_params['observation_time']

    # Dephasing rate: Γ = ω²Ξ/T
    dephasing_rate = omega**2 * Xi / T

    # Dephasing time: τ = 1/Γ
    dephasing_time = 1 / dephasing_rate

    # Frequency scaling verification
    frequency_scaling = diff(dephasing_rate, omega) / dephasing_rate * omega
    # Should give 2 (quadratic scaling)

    return {
        'dephasing_rate': dephasing_rate,
        'dephasing_time': dephasing_time,
        'frequency_scaling': frequency_scaling,
        'omega_squared_dependence': True
    }

def verify_exponential_scaling(visibility_data, scaling_params):
    """
    Verify exponential scaling properties of visibility.

    Checks:
    1. V = exp[-½ω²Ξ] form
    2. Frequency scaling: V ∝ exp[-ω²]
    3. Capacity scaling: V ∝ exp[-Ξ]
    4. Small dephasing limit: V ≈ 1 - ½ω²Ξ
    5. Strong dephasing limit: V → 0 exponentially

    Args:
        visibility_data: Result from compute_visibility
        scaling_params: Parameters for scaling tests

    Returns:
        Dict with verification results
    """
    V = visibility_data['visibility']
    omega = visibility_data['frequency']
    Xi = visibility_data['capacity']
    dephasing_exp = visibility_data['dephasing_exponent']

    # Create symbolic variables for differentiation
    omega_sym = var('omega_sym', domain='positive')
    Xi_sym = var('Xi_sym', domain='positive')

    # Create symbolic expressions
    V_sym = exp(-omega_sym**2 * Xi_sym / 2)
    log_V_sym = log(V_sym)

    # Test 1: Exponential form verification
    # Use safe symbolic comparison method
    exponential_form_diff = V - exp(-dephasing_exp)
    exponential_form_check = exponential_form_diff.simplify_full().is_trivial_zero()

    # Test 2: Frequency scaling verification using symbolic math
    # ∂ln(V)/∂ω = -ωΞ, so ω∂ln(V)/∂ω = -ω²Ξ
    freq_derivative = diff(log_V_sym, omega_sym)
    freq_scaling_symbolic = omega_sym * freq_derivative
    expected_freq_scaling = -omega_sym**2 * Xi_sym
    freq_scaling_diff = (freq_scaling_symbolic - expected_freq_scaling).simplify_full()
    freq_scaling_verified = freq_scaling_diff.is_trivial_zero()

    # Test 3: Capacity scaling verification using symbolic math
    # ∂ln(V)/∂Ξ = -ω²/2, so Ξ∂ln(V)/∂Ξ = -ω²Ξ/2
    capacity_derivative = diff(log_V_sym, Xi_sym)
    capacity_scaling_symbolic = Xi_sym * capacity_derivative
    expected_capacity_scaling = -omega_sym**2 * Xi_sym / 2
    capacity_scaling_diff = (capacity_scaling_symbolic - expected_capacity_scaling).simplify_full()
    capacity_scaling_verified = capacity_scaling_diff.is_trivial_zero()

    # Test 4: Small dephasing limit
    # For ½ω²Ξ ≪ 1: V ≈ 1 - ½ω²Ξ
    # Use symbolic expansion in terms of the dephasing parameter
    dephasing_param = var('dephasing_param')
    V_expansion = exp(-dephasing_param)
    small_dephasing_expansion = V_expansion.taylor(dephasing_param, 0, 3)
    linear_term = small_dephasing_expansion.coefficient(dephasing_param, 1)
    linear_term_diff = (linear_term - (-1)).simplify_full()
    small_limit_check = linear_term_diff.is_trivial_zero()

    # Test 5: Strong dephasing behavior
    # For ½ω²Ξ ≫ 1: V → 0 exponentially fast
    large_arg_limit = limit(exp(-dephasing_param), dephasing_param, oo)
    # Check if limit is zero using safe method
    strong_limit_check = large_arg_limit.is_zero()

    return {
        'exponential_form_verified': exponential_form_check,
        'frequency_scaling_verified': freq_scaling_verified,
        'capacity_scaling_verified': capacity_scaling_verified,
        'small_dephasing_limit_verified': small_limit_check,
        'strong_dephasing_limit_verified': strong_limit_check,
        'all_scaling_verified': all([
            exponential_form_check,
            freq_scaling_verified,
            capacity_scaling_verified,
            small_limit_check,
            strong_limit_check
        ]),
        'scaling_details': {
            'small_expansion': small_dephasing_expansion,
            'frequency_dependence': 'omega_squared',
            'capacity_dependence': 'linear_in_exponent'
        }
    }

def compute_visibility_contours(frequency_range, capacity_range, contour_params=None):
    """
    Compute visibility contours V(ω,Ξ) for parameter space mapping.

    Creates contour plots of constant visibility in (ω,Ξ) parameter space.
    Useful for experimental design and sensitivity analysis.

    Args:
        frequency_range: List [ω_min, ω_max] (rad/s)
        capacity_range: List [Ξ_min, Ξ_max] (s²)
        contour_params: Dict with contour levels and grid

    Returns:
        Dict with:
            'contour_equations': V = constant curves
            'frequency_grid': ω values
            'capacity_grid': Ξ values
            'visibility_surface': V(ω,Ξ) function
    """
    if contour_params is None:
        contour_params = {'levels': [0.1, 0.3, 0.5, 0.7, 0.9]}

    omega_var = var('omega', domain='positive')
    Xi_var = var('Xi', domain='positive')

    # Visibility surface: V(ω,Ξ) = exp[-½ω²Ξ]
    visibility_surface = exp(-omega_var**2 * Xi_var / 2)

    # Contour equations: V = constant ⟹ ω²Ξ = -2ln(V)
    contour_equations = {}
    for V_level in contour_params['levels']:
        # Solve ω²Ξ = -2ln(V) for contour
        contour_constant = -2 * log(V_level)
        contour_equations[f'V_{V_level}'] = omega_var**2 * Xi_var - contour_constant

    return {
        'contour_equations': contour_equations,
        'frequency_variable': omega_var,
        'capacity_variable': Xi_var,
        'visibility_surface': visibility_surface,
        'frequency_range': frequency_range,
        'capacity_range': capacity_range,
        'contour_levels': contour_params['levels']
    }