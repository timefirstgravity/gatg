#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Standard Kerr Solution

Implements the standard GR approach: rotating metric ansatz → Einstein equations → Kerr solution
Computes actual symbolic expressions for Kerr spacetime in Boyer-Lindquist coordinates.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *
from core import setup_manifold_and_coordinates, compute_einstein_tensor, extract_einstein_components

def compute_kerr_functions(r, theta, a, M):
    """
    Compute auxiliary functions for Kerr metric

    Args:
        r, theta: Coordinates
        a: Angular momentum parameter
        M: Mass parameter

    Returns:
        tuple: (Sigma, Delta) where Sigma = r² + a²cos²θ, Delta = r² - 2Mr + a²
    """
    Sigma = r**2 + a**2 * cos(theta)**2
    Delta = r**2 - 2*M*r + a**2
    return Sigma, Delta

def compute_standard_kerr_setup():
    """
    Set up Kerr metric components in Boyer-Lindquist coordinates

    Returns:
        dict: Kerr metric setup and components
    """
    M, X, t, r, th, ph = setup_manifold_and_coordinates()
    assume(r > 0)

    # Kerr parameters
    M_mass = var('M', domain='positive')
    a = var('a', domain='real')  # Angular momentum parameter

    # Auxiliary functions
    Sigma, Delta = compute_kerr_functions(r, th, a, M_mass)

    # Kerr metric components
    g_tt = -(1 - 2*M_mass*r/Sigma)
    g_tph = -2*a*M_mass*r*sin(th)**2/Sigma
    g_rr = Sigma/Delta
    g_thth = Sigma
    g_phph = sin(th)**2 * ((r**2 + a**2) + 2*M_mass*a**2*r*sin(th)**2/Sigma)

    return {
        'manifold': M,
        'coordinates': [t, r, th, ph],
        'mass_parameter': M_mass,
        'angular_parameter': a,
        'auxiliary_functions': {'Sigma': Sigma, 'Delta': Delta},
        'metric_components': {
            'g_tt': g_tt,
            'g_tph': g_tph,
            'g_rr': g_rr,
            'g_thth': g_thth,
            'g_phph': g_phph
        }
    }

def construct_kerr_metric(setup_data):
    """
    Construct complete Kerr metric tensor

    Args:
        setup_data: Result from compute_standard_kerr_setup()

    Returns:
        dict: Complete Kerr metric and verification data
    """
    M = setup_data['manifold']
    components = setup_data['metric_components']

    # Create metric tensor
    g = M.metric('g')

    # Set metric components
    g[0,0] = components['g_tt']     # g_tt
    g[0,3] = components['g_tph']    # g_tφ
    g[3,0] = components['g_tph']    # g_φt (symmetric)
    g[1,1] = components['g_rr']     # g_rr
    g[2,2] = components['g_thth']   # g_θθ
    g[3,3] = components['g_phph']   # g_φφ

    return {
        'metric': g,
        'is_vacuum_solution': True,
        'coordinate_system': 'Boyer-Lindquist',
        'rotation_parameter': setup_data['angular_parameter']
    }

def verify_kerr_einstein_equations(setup_data, metric_data):
    """
    Verify Kerr metric satisfies Einstein vacuum equations

    Args:
        setup_data: Kerr setup data
        metric_data: Constructed metric

    Returns:
        dict: Einstein equation verification results
    """
    g = metric_data['metric']

    # Compute Einstein tensor
    G, R_tensor, R_scalar = compute_einstein_tensor(g)
    components = extract_einstein_components(G, simplify=True)

    # Check all components are zero (vacuum solution)
    all_zero = all(comp.simplify_full() == 0 for comp in components.values())

    return {
        'einstein_components': components,
        'vacuum_verified': all_zero,
        'ricci_scalar': R_scalar,
        'einstein_tensor': G
    }

def verify_schwarzschild_limit(setup_data):
    """
    Verify Kerr reduces to Schwarzschild when a = 0

    Args:
        setup_data: Kerr setup data

    Returns:
        dict: Schwarzschild limit verification
    """
    components = setup_data['metric_components']
    a = setup_data['angular_parameter']

    # Take a → 0 limit
    schwarzschild_components = {}
    for key, comp in components.items():
        schwarzschild_components[key] = comp.substitute(a=0).simplify_full()

    # Expected Schwarzschild results
    r = setup_data['coordinates'][1]
    M_mass = setup_data['mass_parameter']
    rs = 2*M_mass

    expected_g_tt = -(1 - rs/r)
    expected_g_rr = 1/(1 - rs/r)

    # Verify components
    g_tt_correct = (schwarzschild_components['g_tt'] - expected_g_tt).simplify_full() == 0
    g_rr_correct = (schwarzschild_components['g_rr'] - expected_g_rr).simplify_full() == 0
    g_tph_zero = schwarzschild_components['g_tph'].simplify_full() == 0

    limit_correct = g_tt_correct and g_rr_correct and g_tph_zero

    return {
        'schwarzschild_components': schwarzschild_components,
        'limit_verification': limit_correct,
        'frame_dragging_vanishes': g_tph_zero,
        'temporal_component_correct': g_tt_correct,
        'radial_component_correct': g_rr_correct
    }

def compute_slow_rotation_expansion(setup_data):
    """
    Compute first-order slow-rotation expansion of Kerr metric

    Expands metric components in powers of a/M to first order.
    This approximation is valid when |a| << M (slow rotation).

    Args:
        setup_data: Kerr setup from compute_standard_kerr_setup()

    Returns:
        dict: First-order expanded metric components and frame-dragging effects
    """
    components = setup_data['metric_components']
    a = setup_data['angular_parameter']
    M_mass = setup_data['mass_parameter']
    r = setup_data['coordinates'][1]
    th = setup_data['coordinates'][2]

    # Expand each component to first order in a
    expanded_components = {}

    for key, comp in components.items():
        # Use Taylor expansion around a = 0 to first order
        try:
            # Try built-in series expansion
            expanded = comp.series(a, 0, 2).truncate()
        except:
            # Manual expansion: f(a) ≈ f(0) + f'(0)·a
            f_0 = comp.subs({a: 0})
            f_prime_0 = diff(comp, a).subs({a: 0})
            expanded = f_0 + f_prime_0 * a

        expanded_components[key] = expanded.simplify_full()

    # Extract specific first-order terms
    # g_tφ ~ -2Ma*sin²θ/r (frame-dragging term)
    try:
        frame_dragging_coeff = expanded_components['g_tph'].coefficient(a, 1)
    except:
        # Extract manually if coefficient method fails
        g_tph_expanded = expanded_components['g_tph']
        frame_dragging_coeff = diff(g_tph_expanded, a).subs({a: 0})

    expected_frame_dragging = -2*M_mass*sin(th)**2/r

    # g_tt ~ -(1 - 2M/r) + O(a²)
    try:
        g_tt_zeroth = expanded_components['g_tt'].coefficient(a, 0)
    except:
        g_tt_zeroth = expanded_components['g_tt'].subs({a: 0})

    expected_g_tt_zeroth = -(1 - 2*M_mass/r)

    # g_rr ~ (1 - 2M/r)⁻¹ + O(a²)
    try:
        g_rr_zeroth = expanded_components['g_rr'].coefficient(a, 0)
    except:
        g_rr_zeroth = expanded_components['g_rr'].subs({a: 0})

    expected_g_rr_zeroth = 1/(1 - 2*M_mass/r)

    # Verify expansions match expected slow-rotation physics
    frame_dragging_correct = (frame_dragging_coeff - expected_frame_dragging).simplify_full() == 0
    g_tt_correct = (g_tt_zeroth - expected_g_tt_zeroth).simplify_full() == 0
    g_rr_correct = (g_rr_zeroth - expected_g_rr_zeroth).simplify_full() == 0

    return {
        'expanded_components': expanded_components,
        'frame_dragging_coefficient': frame_dragging_coeff,
        'expected_frame_dragging': expected_frame_dragging,
        'frame_dragging_verified': frame_dragging_correct,
        'temporal_zeroth_order': g_tt_zeroth,
        'radial_zeroth_order': g_rr_zeroth,
        'temporal_expansion_correct': g_tt_correct,
        'radial_expansion_correct': g_rr_correct,
        'expansion_order': 'first_order_in_a'
    }

def verify_frame_dragging_physics(setup_data, expansion_data):
    """
    Verify that frame-dragging effects match expected physics

    Frame-dragging (Lense-Thirring effect) should:
    1. Be proportional to angular momentum (∝ a)
    2. Have sin²θ angular dependence (maximal at equator)
    3. Fall off as 1/r (gravitomagnetic field)

    Args:
        setup_data: Kerr metric setup
        expansion_data: Results from slow-rotation expansion

    Returns:
        dict: Frame-dragging physics verification
    """
    frame_dragging = expansion_data['frame_dragging_coefficient']
    a = setup_data['angular_parameter']
    M_mass = setup_data['mass_parameter']
    r = setup_data['coordinates'][1]
    th = setup_data['coordinates'][2]

    # Extract dependences
    # Check proportionality to a (angular momentum)
    a_dependence = frame_dragging.coefficient(a, 1) if frame_dragging.has(a) else frame_dragging

    # Check 1/r falloff
    r_power = -1  # Expected r⁻¹ dependence
    has_correct_r_dependence = a_dependence.has(r**r_power) or a_dependence.has(1/r)

    # Check sin²θ angular dependence
    has_sin2_theta = a_dependence.has(sin(th)**2)

    # Expected magnitude: 2Ma/r for Kerr
    expected_magnitude = 2*M_mass/r
    magnitude_match = (abs(a_dependence) - expected_magnitude*sin(th)**2).simplify_full() == 0

    # Physical interpretation
    lense_thirring_frequency = 2*M_mass*a/(r**3)  # Frame-dragging frequency Ω_LT

    physics_verified = has_correct_r_dependence and has_sin2_theta and magnitude_match

    return {
        'frame_dragging_term': frame_dragging,
        'r_dependence_correct': has_correct_r_dependence,
        'angular_dependence_correct': has_sin2_theta,
        'magnitude_correct': magnitude_match,
        'lense_thirring_frequency': lense_thirring_frequency,
        'physics_verified': physics_verified,
        'physical_meaning': 'Frame-dragging causes test particles to precess with Ω_LT'
    }