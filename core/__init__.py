#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Core Framework
Common functions for all GATG derivations (flux law, Schwarzschild, cosmology, etc.)
"""

from .manifolds import setup_manifold_and_coordinates, setup_static_manifold_and_coordinates
from .metrics import construct_metric, define_temporal_potential, define_static_temporal_potential, compute_3metric_components
from .tensors import compute_einstein_tensor, compute_ricci_tensor, extract_einstein_components, verify_vacuum_solution, compute_einstein_tensor_component
from .constraints import compute_hamiltonian_constraint, compute_general_lambda_momentum_constraint
from .adm import compute_extrinsic_curvature, compute_mixed_extrinsic_curvature, compute_Y_tensor, compute_covariant_divergence
from .symbolic_equations import equations, get_string, get_symbolic
from .differential_operators import dalembertian, spatial_laplacian, gradient_3d, divergence_3d, time_derivative

__all__ = [
    'setup_manifold_and_coordinates',
    'setup_static_manifold_and_coordinates',
    'construct_metric',
    'define_temporal_potential',
    'define_static_temporal_potential',
    'compute_3metric_components',
    'compute_einstein_tensor',
    'compute_ricci_tensor',
    'extract_einstein_components',
    'verify_vacuum_solution',
    'compute_einstein_tensor_component',
    'compute_hamiltonian_constraint',
    'compute_general_lambda_momentum_constraint',
    'compute_extrinsic_curvature',
    'compute_mixed_extrinsic_curvature',
    'compute_Y_tensor',
    'compute_covariant_divergence',
    'equations',
    'get_string',
    'get_symbolic',
    'dalembertian',
    'spatial_laplacian',
    'gradient_3d',
    'divergence_3d',
    'time_derivative'
]