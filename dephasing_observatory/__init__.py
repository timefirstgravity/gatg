#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Dephasing Observatory Module

Gravitational decoherence and capacity calculations.
Implements the dephasing predictions from GATG Quantum Origin paper.

This module computes the capacity Ξ, visibility V, and linewidth Δν
for gravitational dephasing effects in clock networks. All computations
use symbolic SageMath for mathematical precision.
"""

from .capacity import (
    compute_capacity_xi,
    compute_windowed_capacity,
    verify_capacity_properties
)

from .visibility import (
    compute_visibility,
    compute_linewidth,
    compute_dephasing_rate,
    verify_exponential_scaling
)

from .snr_analysis import (
    compute_fisher_matrix,
    compute_parameter_errors,
    build_snr_isolines,
    verify_cramer_rao_bound
)

from .allan_variance import (
    compute_allan_deviation,
    convert_psd_to_allan,
    verify_allan_scaling
)

from .scenarios import (
    leo_ground_scenario,
    ground_ground_scenario,
    intercontinental_scenario,
    verify_physical_constraints
)

__all__ = [
    # Capacity computations
    'compute_capacity_xi',
    'compute_windowed_capacity',
    'verify_capacity_properties',

    # Visibility and dephasing
    'compute_visibility',
    'compute_linewidth',
    'compute_dephasing_rate',
    'verify_exponential_scaling',

    # SNR and parameter estimation
    'compute_fisher_matrix',
    'compute_parameter_errors',
    'build_snr_isolines',
    'verify_cramer_rao_bound',

    # Allan variance bridge
    'compute_allan_deviation',
    'convert_psd_to_allan',
    'verify_allan_scaling',

    # Physical scenarios
    'leo_ground_scenario',
    'ground_ground_scenario',
    'intercontinental_scenario',
    'verify_physical_constraints'
]