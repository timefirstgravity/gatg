#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Standard Gravitoelectromagnetism (GEM)

Implements the standard GR approach: weak field limit → gravitoelectric/gravitomagnetic fields
Computes actual symbolic expressions for GEM field equations and forces.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

def compute_weak_field_metric():
    """
    Set up weak field metric perturbations for GEM derivation

    Returns:
        dict: Weak field metric setup
    """
    # Coordinates
    t, x, y, z = var('t x y z')
    coordinates = [t, x, y, z]

    # Weak field metric: g_μν = η_μν + h_μν where |h_μν| << 1
    # Gravitoelectric potential: Φ (analogous to electric potential)
    # Gravitomagnetic vector potential: A_i (analogous to magnetic vector potential)

    Phi = function('Phi')(t, x, y, z)  # Gravitoelectric potential
    A_x = function('A_x')(t, x, y, z)  # Gravitomagnetic vector potential
    A_y = function('A_y')(t, x, y, z)
    A_z = function('A_z')(t, x, y, z)

    # Metric perturbations in terms of potentials
    # h_00 = -2Φ/c²
    # h_0i = -4A_i/c
    # h_ij = 2Φδ_ij/c² (for static case)

    c = var('c')  # Speed of light

    h_tt = -2*Phi/c**2
    h_tx = -4*A_x/c
    h_ty = -4*A_y/c
    h_tz = -4*A_z/c
    h_xx = 2*Phi/c**2
    h_yy = 2*Phi/c**2
    h_zz = 2*Phi/c**2

    return {
        'coordinates': coordinates,
        'gravitoelectric_potential': Phi,
        'gravitomagnetic_vector': [A_x, A_y, A_z],
        'metric_perturbations': {
            'h_tt': h_tt,
            'h_tx': h_tx,
            'h_ty': h_ty,
            'h_tz': h_tz,
            'h_xx': h_xx,
            'h_yy': h_yy,
            'h_zz': h_zz
        },
        'speed_of_light': c
    }

def derive_gem_fields(weak_field_data):
    """
    Derive gravitoelectric and gravitomagnetic field definitions

    Args:
        weak_field_data: Weak field metric setup

    Returns:
        dict: GEM field definitions
    """
    t, x, y, z = weak_field_data['coordinates']
    Phi = weak_field_data['gravitoelectric_potential']
    A_x, A_y, A_z = weak_field_data['gravitomagnetic_vector']
    c = weak_field_data['speed_of_light']

    # Gravitoelectric field: E_i = -∇_i Φ - (1/c) ∂_t A_i
    E_x = -diff(Phi, x) - (1/c) * diff(A_x, t)
    E_y = -diff(Phi, y) - (1/c) * diff(A_y, t)
    E_z = -diff(Phi, z) - (1/c) * diff(A_z, t)

    # Gravitomagnetic field: B_i = (∇ × A)_i
    B_x = diff(A_z, y) - diff(A_y, z)
    B_y = diff(A_x, z) - diff(A_z, x)
    B_z = diff(A_y, x) - diff(A_x, y)

    return {
        'gravitoelectric_field': [E_x, E_y, E_z],
        'gravitomagnetic_field': [B_x, B_y, B_z],
        'field_definitions': {
            'E_field': 'E = -∇Φ - (1/c)∂A/∂t',
            'B_field': 'B = ∇ × A'
        }
    }

def derive_gem_field_equations(gem_fields, weak_field_data):
    """
    Derive GEM field equations from linearized Einstein equations

    Args:
        gem_fields: GEM field definitions
        weak_field_data: Weak field setup

    Returns:
        dict: GEM field equations (Maxwell-like)
    """
    t, x, y, z = weak_field_data['coordinates']
    E_x, E_y, E_z = gem_fields['gravitoelectric_field']
    B_x, B_y, B_z = gem_fields['gravitomagnetic_field']
    c = weak_field_data['speed_of_light']

    # Energy-momentum source terms
    rho = var('rho')      # Energy density
    j_x = var('j_x')      # Energy flux density
    j_y = var('j_y')
    j_z = var('j_z')
    G = var('G')          # Gravitational constant

    # GEM field equations (analogous to Maxwell equations)
    # 1. Gauss law for gravitoelectricity: ∇·E = -4πGρ/c²
    gauss_law = diff(E_x, x) + diff(E_y, y) + diff(E_z, z) + 4*pi*G*rho/c**2

    # 2. No gravitomagnetic monopoles: ∇·B = 0
    no_monopoles = diff(B_x, x) + diff(B_y, y) + diff(B_z, z)

    # 3. Faraday law for gravitomagnetism: ∇×E = -(1/c)∂B/∂t
    faraday_x = diff(E_z, y) - diff(E_y, z) + (1/c) * diff(B_x, t)
    faraday_y = diff(E_x, z) - diff(E_z, x) + (1/c) * diff(B_y, t)
    faraday_z = diff(E_y, x) - diff(E_x, y) + (1/c) * diff(B_z, t)

    # 4. Ampère law for gravitomagnetism: ∇×B = -(4πG/c²)j - (1/c)∂E/∂t
    ampere_x = diff(B_z, y) - diff(B_y, z) + (4*pi*G/c**2)*j_x + (1/c) * diff(E_x, t)
    ampere_y = diff(B_x, z) - diff(B_z, x) + (4*pi*G/c**2)*j_y + (1/c) * diff(E_y, t)
    ampere_z = diff(B_y, x) - diff(B_x, y) + (4*pi*G/c**2)*j_z + (1/c) * diff(E_z, t)

    return {
        'gauss_law': gauss_law,
        'no_monopoles': no_monopoles,
        'faraday_laws': [faraday_x, faraday_y, faraday_z],
        'ampere_laws': [ampere_x, ampere_y, ampere_z],
        'source_terms': {
            'energy_density': rho,
            'energy_flux': [j_x, j_y, j_z]
        },
        'field_equations_form': 'Maxwell-like equations for gravity'
    }

def compute_gem_forces(gem_fields):
    """
    Compute gravitational forces on test masses in GEM formulation

    Args:
        gem_fields: GEM field definitions

    Returns:
        dict: Gravitational force expressions
    """
    E_x, E_y, E_z = gem_fields['gravitoelectric_field']
    B_x, B_y, B_z = gem_fields['gravitomagnetic_field']

    # Test mass variables
    m = var('m')          # Test mass
    v_x = var('v_x')      # Test mass velocity
    v_y = var('v_y')
    v_z = var('v_z')
    c = var('c')

    # Gravitational force: F = m[E + v×B] (analogous to Lorentz force)
    # Gravitoelectric force: F_E = mE
    F_E_x = m * E_x
    F_E_y = m * E_y
    F_E_z = m * E_z

    # Gravitomagnetic force: F_B = m(v × B)
    F_B_x = m * (v_y * B_z - v_z * B_y)
    F_B_y = m * (v_z * B_x - v_x * B_z)
    F_B_z = m * (v_x * B_y - v_y * B_x)

    # Total force
    F_total_x = F_E_x + F_B_x
    F_total_y = F_E_y + F_B_y
    F_total_z = F_E_z + F_B_z

    return {
        'gravitoelectric_force': [F_E_x, F_E_y, F_E_z],
        'gravitomagnetic_force': [F_B_x, F_B_y, F_B_z],
        'total_force': [F_total_x, F_total_y, F_total_z],
        'force_law': 'F = m(E + v × B)',
        'test_mass': m,
        'test_velocity': [v_x, v_y, v_z]
    }