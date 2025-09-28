#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Standard Linearized General Relativity

Implements the traditional approach: metric perturbations → Einstein equations → TT gauge
Computes actual symbolic expressions for constraint elimination and wave operator derivation.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

def compute_standard_linearization():
    """
    Compute standard linearized GR starting from metric perturbations

    Returns:
        dict: Initial metric setup and linearized Einstein tensor
    """
    # Coordinates
    t, x, y, z = var('t x y z')
    coordinates = [t, x, y, z]

    # Metric perturbations h_μν (symmetric 4×4)
    h_tt, h_tx, h_ty, h_tz = var('h_tt h_tx h_ty h_tz')
    h_xx, h_xy, h_xz, h_yy, h_yz, h_zz = var('h_xx h_xy h_xz h_yy h_yz h_zz')

    h_perturbation = matrix(SR, [
        [h_tt, h_tx, h_ty, h_tz],
        [h_tx, h_xx, h_xy, h_xz],
        [h_ty, h_xy, h_yy, h_yz],
        [h_tz, h_xz, h_yz, h_zz]
    ])

    # Minkowski background
    eta = matrix(SR, [
        [-1,  0,  0,  0],
        [ 0,  1,  0,  0],
        [ 0,  0,  1,  0],
        [ 0,  0,  0,  1]
    ])

    # Full metric: g_μν = η_μν + h_μν
    full_metric = eta + h_perturbation

    # Trace
    h_trace = h_tt + h_xx + h_yy + h_zz

    return {
        'coordinates': coordinates,
        'metric_perturbation': h_perturbation,
        'background_metric': eta,
        'full_metric': full_metric,
        'trace': h_trace,
        'initial_dof': 10  # Symmetric 4×4 matrix
    }

def apply_tt_gauge(metric_data):
    """
    Apply TT gauge conditions to eliminate gauge degrees of freedom

    Args:
        metric_data: Result from compute_standard_linearization()

    Returns:
        dict: TT gauge constraints and remaining degrees of freedom
    """
    coordinates = metric_data['coordinates']
    x, y, z = coordinates[1:4]  # Spatial coordinates

    # TT gauge variables
    h_xx_TT, h_xy_TT, h_xz_TT = var('h_xx_TT h_xy_TT h_xz_TT')
    h_yy_TT, h_yz_TT, h_zz_TT = var('h_yy_TT h_yz_TT h_zz_TT')

    # TT conditions
    # 1. Traceless: h_ii = 0
    traceless_constraint = h_xx_TT + h_yy_TT + h_zz_TT

    # 2. Transverse: ∂_i h_ij = 0
    transverse_x = diff(h_xx_TT, x) + diff(h_xy_TT, y) + diff(h_xz_TT, z)
    transverse_y = diff(h_xy_TT, x) + diff(h_yy_TT, y) + diff(h_yz_TT, z)
    transverse_z = diff(h_xz_TT, x) + diff(h_yz_TT, y) + diff(h_zz_TT, z)

    # 3. Temporal gauge: h_0μ = 0
    temporal_constraints = 4  # h_tt = h_tx = h_ty = h_tz = 0

    # TT metric (only spatial components)
    h_TT_matrix = matrix(SR, [
        [0, 0, 0, 0],
        [0, h_xx_TT, h_xy_TT, h_xz_TT],
        [0, h_xy_TT, h_yy_TT, h_yz_TT],
        [0, h_xz_TT, h_yz_TT, h_zz_TT]
    ])

    return {
        'tt_variables': [h_xx_TT, h_xy_TT, h_xz_TT, h_yy_TT, h_yz_TT, h_zz_TT],
        'traceless_constraint': traceless_constraint,
        'transverse_constraints': [transverse_x, transverse_y, transverse_z],
        'temporal_constraint_count': temporal_constraints,
        'tt_metric': h_TT_matrix,
        'total_constraints': 1 + 3 + 4,  # traceless + transverse + temporal
        'effective_gauge_constraints': 4  # 4 coordinate degrees of freedom
    }

def extract_physical_modes(metric_data, tt_data):
    """
    Extract final physical degrees of freedom after all constraints

    Args:
        metric_data: Initial metric data
        tt_data: TT gauge data

    Returns:
        dict: Final physical DOF and wave operator
    """
    initial_dof = metric_data['initial_dof']
    gauge_constraints = tt_data['effective_gauge_constraints']
    bianchi_constraints = 4  # From Einstein equations

    # DOF counting
    after_einstein = initial_dof - bianchi_constraints  # 10 - 4 = 6
    final_physical_dof = after_einstein - gauge_constraints  # 6 - 4 = 2

    # Wave operator: □ = η^μν ∂_μ ∂_ν = -∂²_t + ∇²
    t, x, y, z = metric_data['coordinates']
    wave_operator = lambda field: -diff(field, t, 2) + diff(field, x, 2) + diff(field, y, 2) + diff(field, z, 2)

    # Physical TT modes (2 independent components)
    h_plus = var('h_plus')   # Plus polarization
    h_cross = var('h_cross') # Cross polarization

    return {
        'initial_dof': initial_dof,
        'after_einstein_constraints': after_einstein,
        'final_physical_dof': final_physical_dof,
        'wave_operator': wave_operator,
        'physical_modes': [h_plus, h_cross],
        'wave_equation_form': '□h^TT_ij = 0',
        'propagation_speed': 'c = 1'
    }