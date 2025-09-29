#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Quantum Origin Module

Implements the fixed-point emergence of classical spacetime from quantum redundancy.
Based on GATG Quantum Origin paper: Φ = T[Φ] ≡ C ω² C[Φ].

This module demonstrates how the lapse function Φ emerges as a fixed point
of the operator T constructed from the CTP (Closed Time Path) kernel C.
The contraction property |T| < 1 ensures well-posed emergence.
"""

from .ctp_kernel import (
    build_local_kernel,
    build_quasilocal_kernel,
    build_screened_kernel,
    compute_kernel_spectrum,
    verify_kernel_properties
)

from .fixed_point import (
    build_fixed_point_operator,
    check_contraction_property,
    iterate_fixed_point,
    compute_spectral_radius,
    verify_banach_conditions
)

from .screening_length import (
    extract_screening_length,
    compute_correlation_function,
    analyze_greens_function_tail,
    verify_screening_bounds
)

from .capacity_functional import (
    compute_capacity_from_kernel,
    build_capacity_functional,
    compute_ctp_noise_spectrum,
    verify_fdt_relation
)

from .physical_parameters import (
    get_benchmark_parameters,
    solve_contraction_constraints,
    benchmark_fixed_point_viability
)

__all__ = [
    # CTP kernel functions
    'build_local_kernel',
    'build_quasilocal_kernel',
    'build_screened_kernel',
    'compute_kernel_spectrum',
    'verify_kernel_properties',

    # Fixed-point operations
    'build_fixed_point_operator',
    'check_contraction_property',
    'iterate_fixed_point',
    'compute_spectral_radius',
    'verify_banach_conditions',

    # Screening length analysis
    'extract_screening_length',
    'compute_correlation_function',
    'analyze_greens_function_tail',
    'verify_screening_bounds',

    # Capacity functional
    'compute_capacity_from_kernel',
    'build_capacity_functional',
    'compute_ctp_noise_spectrum',
    'verify_fdt_relation',

    # Physical parameters
    'get_benchmark_parameters',
    'solve_contraction_constraints',
    'benchmark_fixed_point_viability'
]