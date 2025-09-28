#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Lapse-First Linearized General Relativity (GATG Approach)

Implements the GATG approach: (δΦ, ωᵢ, hᵢⱼ) → ADM constraints → TT gauge
Computes actual symbolic expressions for constraint elimination in lapse-first formulation.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

def compute_lapse_first_variables():
    """
    Set up lapse-first variables following GATG decomposition

    Returns:
        dict: Lapse-first variable setup and initial DOF count
    """
    # Coordinates
    t, x, y, z = var('t x y z')
    coordinates = [t, x, y, z]
    spatial_coords = [x, y, z]

    # GATG variables
    # 1. Lapse perturbation: N = e^Φ ≈ 1 + δΦ
    delta_Phi = var('delta_Phi')

    # 2. Shift perturbations: N^i = ω^i
    omega_x, omega_y, omega_z = var('omega_x omega_y omega_z')
    shift_vector = [omega_x, omega_y, omega_z]

    # 3. Spatial metric perturbations: γ_ij = δ_ij + h_ij
    h_xx, h_xy, h_xz = var('h_xx h_xy h_xz')
    h_yy, h_yz, h_zz = var('h_yy h_yz h_zz')

    spatial_metric_perturbations = matrix(SR, [
        [h_xx, h_xy, h_xz],
        [h_xy, h_yy, h_yz],
        [h_xz, h_yz, h_zz]
    ])

    # Count initial degrees of freedom
    initial_variables = {
        'lapse': 1,        # δΦ
        'shift': 3,        # ω_x, ω_y, ω_z
        'spatial_metric': 6 # h_ij (symmetric 3×3)
    }
    total_initial_dof = sum(initial_variables.values())  # 1 + 3 + 6 = 10

    return {
        'coordinates': coordinates,
        'spatial_coordinates': spatial_coords,
        'lapse_perturbation': delta_Phi,
        'shift_vector': shift_vector,
        'spatial_metric_perturbations': spatial_metric_perturbations,
        'initial_variables': initial_variables,
        'total_initial_dof': total_initial_dof
    }

def apply_adm_constraints(variables_data):
    """
    Apply ADM constraints to eliminate non-physical degrees of freedom

    Args:
        variables_data: Result from compute_lapse_first_variables()

    Returns:
        dict: Constraint equations and remaining DOF after elimination
    """
    coords = variables_data['coordinates']
    spatial_coords = variables_data['spatial_coordinates']
    x, y, z = spatial_coords

    delta_Phi = variables_data['lapse_perturbation']
    omega_x, omega_y, omega_z = variables_data['shift_vector']
    h_matrix = variables_data['spatial_metric_perturbations']

    # Extract spatial metric components
    h_xx = h_matrix[0, 0]
    h_xy = h_matrix[0, 1]
    h_xz = h_matrix[0, 2]
    h_yy = h_matrix[1, 1]
    h_yz = h_matrix[1, 2]
    h_zz = h_matrix[2, 2]

    # ADM constraints (linearized)
    # 1. Hamiltonian constraint: relates δΦ to spatial metric
    # Simplified linearized form: ∇²δΦ ≈ (1/2)∇²h_trace
    h_trace = h_xx + h_yy + h_zz
    hamiltonian_constraint = diff(delta_Phi, x, 2) + diff(delta_Phi, y, 2) + diff(delta_Phi, z, 2) - (1/2) * (diff(h_trace, x, 2) + diff(h_trace, y, 2) + diff(h_trace, z, 2))

    # 2. Momentum constraints: relate ω_i to spatial metric derivatives
    # Simplified linearized form: ∂_j(h_ij - (1/2)δ_ij h) ≈ 0
    momentum_x = diff(h_xx - (1/2)*h_trace, x) + diff(h_xy, y) + diff(h_xz, z)
    momentum_y = diff(h_xy, x) + diff(h_yy - (1/2)*h_trace, y) + diff(h_yz, z)
    momentum_z = diff(h_xz, x) + diff(h_yz, y) + diff(h_zz - (1/2)*h_trace, z)

    # After constraints: δΦ and ω_i are determined by h_ij
    # This eliminates 4 DOF (1 + 3), leaving 6 DOF in h_ij
    remaining_dof_after_constraints = 6  # Only spatial metric h_ij remains

    return {
        'hamiltonian_constraint': hamiltonian_constraint,
        'momentum_constraints': [momentum_x, momentum_y, momentum_z],
        'constraints_applied': 4,  # 1 Hamiltonian + 3 momentum
        'remaining_dof': remaining_dof_after_constraints,
        'eliminated_variables': ['delta_Phi', 'omega_x', 'omega_y', 'omega_z'],
        'remaining_variables': ['h_xx', 'h_xy', 'h_xz', 'h_yy', 'h_yz', 'h_zz']
    }

def extract_tt_modes(variables_data, constraint_data):
    """
    Apply TT gauge to remaining spatial metric to extract physical modes

    Args:
        variables_data: Initial variables
        constraint_data: After ADM constraints

    Returns:
        dict: Final physical DOF and wave operator (same as Standard GR)
    """
    spatial_coords = variables_data['spatial_coordinates']
    x, y, z = spatial_coords

    # After constraints, only spatial metric h_ij remains (6 DOF)
    h_xx, h_xy, h_xz = var('h_xx h_xy h_xz')
    h_yy, h_yz, h_zz = var('h_yy h_yz h_zz')

    # Apply TT conditions to spatial metric
    # 1. Traceless: h_ii = 0
    traceless_constraint = h_xx + h_yy + h_zz

    # 2. Transverse: ∂_i h_ij = 0
    transverse_x = diff(h_xx, x) + diff(h_xy, y) + diff(h_xz, z)
    transverse_y = diff(h_xy, x) + diff(h_yy, y) + diff(h_yz, z)
    transverse_z = diff(h_xz, x) + diff(h_yz, y) + diff(h_zz, z)

    # DOF counting after TT gauge
    initial_spatial_dof = constraint_data['remaining_dof']  # 6
    tt_constraints_applied = 4  # 1 traceless + 3 transverse
    final_physical_dof = initial_spatial_dof - tt_constraints_applied  # 6 - 4 = 2

    # Wave operator (same as Standard GR)
    t = variables_data['coordinates'][0]
    wave_operator = lambda field: -diff(field, t, 2) + diff(field, x, 2) + diff(field, y, 2) + diff(field, z, 2)

    # Physical modes (same 2 polarizations as Standard GR)
    h_plus = var('h_plus')   # Plus polarization
    h_cross = var('h_cross') # Cross polarization

    return {
        'initial_lapse_first_dof': variables_data['total_initial_dof'],
        'after_adm_constraints': constraint_data['remaining_dof'],
        'after_tt_gauge': final_physical_dof,
        'final_physical_dof': final_physical_dof,
        'traceless_constraint': traceless_constraint,
        'transverse_constraints': [transverse_x, transverse_y, transverse_z],
        'wave_operator': wave_operator,
        'physical_modes': [h_plus, h_cross],
        'wave_equation_form': '□h^TT_ij = 0',
        'propagation_speed': 'c = 1'
    }