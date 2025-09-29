#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Quantum Witness Module

Clock network quantum signatures and commutator witness computations.
Implements the quantum detection protocols from GATG Paper V.

This module enables detection of quantum signatures in gravitational
dephasing through time-reversal asymmetry and commutator witnesses.
All computations use symbolic SageMath for mathematical rigor.
"""

from .spectral_core import (
    build_psd_lapse_fluctuations,
    build_antisymmetric_spectrum,
    compute_total_filter,
    verify_spectral_properties
)

from .filters import (
    ramsey_filter,
    hahn_echo_filter,
    cpmg_filter,
    baseline_response_function,
    verify_filter_properties
)

from .commutator_witness import (
    compute_q1_witness,
    compute_visibility_difference,
    extract_quantum_signature,
    verify_classical_limit
)

from .cumulants import (
    compute_second_cumulant,
    compute_third_cumulant,
    compute_fourth_cumulant,
    compute_snr_lines
)

__all__ = [
    # Spectral analysis
    'build_psd_lapse_fluctuations',
    'build_antisymmetric_spectrum',
    'compute_total_filter',
    'verify_spectral_properties',

    # Clock sequence filters
    'ramsey_filter',
    'hahn_echo_filter',
    'cpmg_filter',
    'baseline_response_function',
    'verify_filter_properties',

    # Quantum witness computations
    'compute_q1_witness',
    'compute_visibility_difference',
    'extract_quantum_signature',
    'verify_classical_limit',

    # Cumulant analysis
    'compute_second_cumulant',
    'compute_third_cumulant',
    'compute_fourth_cumulant',
    'compute_snr_lines'
]