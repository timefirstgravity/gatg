#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Lapse-First Schwarzschild Solution

Implements the GATG approach: start with temporal potential Φ → construct metric → verify equivalence
Computes actual symbolic expressions for lapse-first formulation of Schwarzschild spacetime.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *
from core import setup_manifold_and_coordinates

def compute_lapse_first_variables():
    """
    Set up lapse-first variables for Schwarzschild spacetime

    Returns:
        dict: Temporal potential Φ and coordinate setup
    """
    M, X, t, r, th, ph = setup_manifold_and_coordinates()
    assume(r > 0)

    # Define symbolic parameters
    rs = var('r_s', domain='positive')

    # Temporal potential Φ where N = e^Φ
    # For Schwarzschild: A(r) = 1 - r_s/r = e^(2Φ)
    # Therefore: Φ = (1/2) ln(1 - r_s/r)
    A_function = 1 - rs/r
    Phi = (1/2) * ln(A_function)

    # Lapse function N = e^Φ
    N = exp(Phi)

    return {
        'manifold': M,
        'coordinates': [t, r, th, ph],
        'schwarzschild_radius': rs,
        'temporal_potential': Phi,
        'lapse_function': N,
        'coordinate_r': r,
        'coordinate_th': th
    }

def construct_spatial_metric(lapse_data):
    """
    Construct spatial 3-metric γᵢⱼ for Schwarzschild spacetime

    Args:
        lapse_data: Result from compute_lapse_first_variables()

    Returns:
        dict: Spatial metric components
    """
    r = lapse_data['coordinate_r']
    th = lapse_data['coordinate_th']
    rs = lapse_data['schwarzschild_radius']

    # Spatial metric components
    # γ_rr = (1 - r_s/r)^(-1)
    # γ_θθ = r²
    # γ_φφ = r² sin²θ
    gamma_rr = 1 / (1 - rs/r)
    gamma_thth = r**2
    gamma_phph = r**2 * sin(th)**2

    return {
        'gamma_rr': gamma_rr,
        'gamma_thth': gamma_thth,
        'gamma_phph': gamma_phph,
        'spatial_dimension': 3
    }

def construct_lapse_first_metric(lapse_data, spatial_data):
    """
    Construct complete 4D metric from lapse function and spatial metric

    Args:
        lapse_data: Lapse function data
        spatial_data: Spatial metric data

    Returns:
        dict: Complete Schwarzschild metric in lapse-first form
    """
    M = lapse_data['manifold']
    N = lapse_data['lapse_function']

    # Create metric tensor
    g = M.metric('g')

    # Lapse-first form: ds² = -N² dt² + γᵢⱼ dx^i dx^j
    g[0,0] = -N**2  # -N² dt²
    g[1,1] = spatial_data['gamma_rr']     # γ_rr dr²
    g[2,2] = spatial_data['gamma_thth']   # γ_θθ dθ²
    g[3,3] = spatial_data['gamma_phph']   # γ_φφ dφ²

    return {
        'metric': g,
        'temporal_component': -N**2,
        'spatial_components': spatial_data,
        'metric_form': 'lapse_first'
    }

def verify_lapse_equivalence(lapse_data, metric_data):
    """
    Verify lapse-first construction gives correct Schwarzschild metric

    Args:
        lapse_data: Lapse function data
        metric_data: Constructed metric

    Returns:
        dict: Verification that ds² = -(1-r_s/r) dt² + ... is recovered
    """
    r = lapse_data['coordinate_r']
    rs = lapse_data['schwarzschild_radius']
    N = lapse_data['lapse_function']

    # Expected Schwarzschild g_tt component
    expected_g_tt = -(1 - rs/r)

    # Computed g_tt from lapse-first: -N²
    computed_g_tt = -N**2

    # Simplify and verify equivalence
    difference = (computed_g_tt - expected_g_tt).simplify_full()
    equivalence_verified = (difference == 0)

    return {
        'expected_g_tt': expected_g_tt,
        'computed_g_tt': computed_g_tt,
        'difference': difference,
        'equivalence_verified': equivalence_verified,
        'lapse_construction_correct': equivalence_verified
    }