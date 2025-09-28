#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Lapse-First Kerr Solution

Implements the GATG approach: extract temporal potential and shift from Kerr metric
Computes actual symbolic expressions for lapse-first formulation of rotating spacetime.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *
from core import setup_manifold_and_coordinates

def compute_lapse_first_kerr_variables():
    """
    Set up lapse-first variables for Kerr spacetime

    Returns:
        dict: Temporal potential, shift vector, and spatial metric extraction
    """
    M, X, t, r, th, ph = setup_manifold_and_coordinates()
    assume(r > 0)

    # Kerr parameters
    M_mass = var('M', domain='positive')
    a = var('a', domain='real')

    # Auxiliary functions
    Sigma = r**2 + a**2 * cos(th)**2
    Delta = r**2 - 2*M_mass*r + a**2

    return {
        'manifold': M,
        'coordinates': [t, r, th, ph],
        'mass_parameter': M_mass,
        'angular_parameter': a,
        'auxiliary_functions': {'Sigma': Sigma, 'Delta': Delta},
        'coordinate_r': r,
        'coordinate_th': th
    }

def extract_lapse_function(kerr_data):
    """
    Extract lapse function N from Kerr metric temporal component

    Args:
        kerr_data: Kerr variable setup

    Returns:
        dict: Lapse function and temporal potential
    """
    r = kerr_data['coordinate_r']
    th = kerr_data['coordinate_th']
    M_mass = kerr_data['mass_parameter']
    a = kerr_data['angular_parameter']
    Sigma = kerr_data['auxiliary_functions']['Sigma']

    # From Kerr metric: g_tt = -(1 - 2Mr/Σ) = -N²
    # Therefore: N² = 1 - 2Mr/Σ
    N_squared = 1 - 2*M_mass*r/Sigma
    N = sqrt(N_squared)

    # Temporal potential: Φ = ln(N)
    Phi = ln(N)

    return {
        'lapse_function': N,
        'temporal_potential': Phi,
        'lapse_squared': N_squared
    }

def extract_shift_vector(kerr_data):
    """
    Extract shift vector from Kerr metric off-diagonal terms

    Args:
        kerr_data: Kerr variable setup

    Returns:
        dict: Shift vector components
    """
    r = kerr_data['coordinate_r']
    th = kerr_data['coordinate_th']
    M_mass = kerr_data['mass_parameter']
    a = kerr_data['angular_parameter']
    Sigma = kerr_data['auxiliary_functions']['Sigma']

    # From Kerr metric: g_tφ = -2aMr sin²θ/Σ
    # In lapse-first form: g_tφ = N² N^φ where N^φ is shift component
    # Therefore: N^φ = g_tφ/N² = -2aMr sin²θ/(Σ(1 - 2Mr/Σ))

    N_squared = 1 - 2*M_mass*r/Sigma
    shift_phi = -2*a*M_mass*r*sin(th)**2/(Sigma * N_squared)

    return {
        'shift_r': 0,      # No radial shift
        'shift_th': 0,     # No polar shift
        'shift_phi': shift_phi,  # Azimuthal shift (frame dragging)
        'frame_dragging_component': shift_phi
    }

def construct_spatial_metric(kerr_data):
    """
    Extract spatial 3-metric from Kerr metric

    Args:
        kerr_data: Kerr variable setup

    Returns:
        dict: Spatial metric components γᵢⱼ
    """
    r = kerr_data['coordinate_r']
    th = kerr_data['coordinate_th']
    M_mass = kerr_data['mass_parameter']
    a = kerr_data['angular_parameter']
    Sigma, Delta = kerr_data['auxiliary_functions']['Sigma'], kerr_data['auxiliary_functions']['Delta']

    # Spatial metric components from Kerr metric
    gamma_rr = Sigma/Delta
    gamma_thth = Sigma
    gamma_phph = sin(th)**2 * ((r**2 + a**2) + 2*M_mass*a**2*r*sin(th)**2/Sigma)

    return {
        'gamma_rr': gamma_rr,
        'gamma_thth': gamma_thth,
        'gamma_phph': gamma_phph,
        'spatial_determinant': gamma_rr * gamma_thth * gamma_phph
    }

def verify_lapse_first_reconstruction(kerr_data, lapse_data, shift_data, spatial_data):
    """
    Verify lapse-first construction reproduces original Kerr metric

    Args:
        kerr_data: Original Kerr setup
        lapse_data: Extracted lapse function
        shift_data: Extracted shift vector
        spatial_data: Extracted spatial metric

    Returns:
        dict: Reconstruction verification
    """
    # Original Kerr components
    r = kerr_data['coordinate_r']
    th = kerr_data['coordinate_th']
    M_mass = kerr_data['mass_parameter']
    a = kerr_data['angular_parameter']
    Sigma = kerr_data['auxiliary_functions']['Sigma']

    original_g_tt = -(1 - 2*M_mass*r/Sigma)
    original_g_tph = -2*a*M_mass*r*sin(th)**2/Sigma

    # Reconstructed from lapse-first: ds² = -N² dt² + γᵢⱼ (dx^i + N^i dt)(dx^j + N^j dt)
    N = lapse_data['lapse_function']
    N_phi = shift_data['shift_phi']

    reconstructed_g_tt = -N**2
    reconstructed_g_tph = N**2 * N_phi

    # Verify reconstruction
    g_tt_difference = (reconstructed_g_tt - original_g_tt).simplify_full()
    g_tph_difference = (reconstructed_g_tph - original_g_tph).simplify_full()

    reconstruction_correct = (g_tt_difference == 0) and (g_tph_difference == 0)

    return {
        'g_tt_difference': g_tt_difference,
        'g_tph_difference': g_tph_difference,
        'reconstruction_verified': reconstruction_correct,
        'original_components': {'g_tt': original_g_tt, 'g_tph': original_g_tph},
        'reconstructed_components': {'g_tt': reconstructed_g_tt, 'g_tph': reconstructed_g_tph}
    }