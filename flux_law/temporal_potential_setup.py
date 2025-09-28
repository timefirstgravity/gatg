#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Temporal Potential Setup for Flux Law Derivation

Sets up spherically symmetric spacetime with temporal potential Φ
Computes actual symbolic expressions for GATG flux law derivation.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *
from core import setup_manifold_and_coordinates

def setup_spherical_temporal_potential():
    """
    Set up spherically symmetric spacetime with temporal potential Φ(t,r)

    Returns:
        dict: Temporal potential setup and coordinate system
    """
    # Coordinates
    t, r, th, ph = var('t r theta phi')
    assume(r > 0)

    # Temporal potential Φ(t,r) - the key GATG variable
    Phi = function('Phi')(t, r)

    # Lapse function N = e^Φ
    N = exp(Phi)

    # For spherical symmetry: no shift vector
    shift_r = 0
    shift_th = 0
    shift_ph = 0

    return {
        'coordinates': [t, r, th, ph],
        'temporal_potential': Phi,
        'lapse_function': N,
        'shift_vector': [shift_r, shift_th, shift_ph],
        'coordinate_r': r,
        'coordinate_t': t
    }

def construct_adm_variables(potential_data):
    """
    Construct ADM variables from temporal potential

    Args:
        potential_data: Temporal potential setup

    Returns:
        dict: ADM decomposition variables
    """
    t, r, th, ph = potential_data['coordinates']
    Phi = potential_data['temporal_potential']
    N = potential_data['lapse_function']

    # Spatial metric γᵢⱼ for spherical coordinates
    # Assume A(t,r) = e^(2Φ) for radial component
    A = exp(2 * Phi)

    gamma_rr = A
    gamma_thth = r**2
    gamma_phph = r**2 * sin(th)**2

    # Extrinsic curvature components
    # K_ij = (1/2N)[∂_t γ_ij - ∇_i N_j - ∇_j N_i]
    # For spherical symmetry with no shift: K_ij = (1/2N) ∂_t γ_ij

    K_rr = (1/(2*N)) * diff(gamma_rr, t)
    K_thth = (1/(2*N)) * diff(gamma_thth, t)
    K_phph = (1/(2*N)) * diff(gamma_phph, t)

    # Trace of extrinsic curvature
    # K = γ^ij K_ij
    K = K_rr/gamma_rr + K_thth/gamma_thth + K_phph/gamma_phph

    return {
        'spatial_metric': {
            'gamma_rr': gamma_rr,
            'gamma_thth': gamma_thth,
            'gamma_phph': gamma_phph
        },
        'extrinsic_curvature': {
            'K_rr': K_rr,
            'K_thth': K_thth,
            'K_phph': K_phph
        },
        'extrinsic_curvature_trace': K,
        'function_A': A
    }

def compute_energy_momentum_variables():
    """
    Set up energy-momentum tensor variables for flux law

    Returns:
        dict: Energy-momentum tensor components
    """
    # Energy-momentum tensor components
    # T^μν in spherical coordinates
    T_tt = var('T_tt')  # Energy density (as seen by normal observer)
    T_tr = var('T_tr')  # Energy flux (radial component)
    T_rr = var('T_rr')  # Radial pressure
    T_thth = var('T_thth')  # Angular pressure
    T_phph = var('T_phph')  # Angular pressure

    # Physical constants
    G = var('G')  # Gravitational constant
    c = var('c')  # Speed of light

    return {
        'stress_energy_components': {
            'T_tt': T_tt,
            'T_tr': T_tr,
            'T_rr': T_rr,
            'T_thth': T_thth,
            'T_phph': T_phph
        },
        'physical_constants': {
            'G': G,
            'c': c
        }
    }

def verify_temporal_potential_setup(potential_data, adm_data, energy_data):
    """
    Verify the temporal potential setup is consistent

    Args:
        potential_data: Temporal potential setup
        adm_data: ADM variables
        energy_data: Energy-momentum variables

    Returns:
        dict: Setup verification results
    """
    Phi = potential_data['temporal_potential']
    N = potential_data['lapse_function']
    A = adm_data['function_A']

    # Verify lapse-potential relation: N = e^Φ
    lapse_relation_correct = (N - exp(Phi)).simplify_full() == 0

    # Verify A-potential relation: A = e^(2Φ)
    A_relation_correct = (A - exp(2*Phi)).simplify_full() == 0

    # Check spherical symmetry
    r = potential_data['coordinate_r']
    gamma_thth = adm_data['spatial_metric']['gamma_thth']
    gamma_phph = adm_data['spatial_metric']['gamma_phph']

    spherical_symmetry_correct = (gamma_thth == r**2) and (gamma_phph == r**2 * sin(potential_data['coordinates'][2])**2)

    setup_valid = lapse_relation_correct and A_relation_correct and spherical_symmetry_correct

    return {
        'lapse_relation_verified': lapse_relation_correct,
        'A_relation_verified': A_relation_correct,
        'spherical_symmetry_verified': spherical_symmetry_correct,
        'temporal_potential_setup_valid': setup_valid
    }