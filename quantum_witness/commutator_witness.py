#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Quantum Witness: Commutator Witness Computation

Implements the Q1 commutator witness protocol from GATG Paper V.
This is the core quantum signature detection mechanism.

The Q1 witness detects quantum non-commutativity through time-reversal
asymmetry in clock dephasing. Classical fields give zero witness,
while quantum fields produce non-zero signatures.
"""

from sage.all import *
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .spectral_core import compute_total_filter
from .filters import verify_filter_properties

def compute_q1_witness(clock_A_filter, clock_B_filter, antisymmetric_psd, frequency_var, integration_params=None):
    """
    Compute Q1 commutator witness.

    From Paper V, Equation Q1:
    Δ(-ln V) = 2 Im ∫ df Ω_A(f) Ω_B*(f) S^-_Φ(f)

    This integral detects quantum non-commutativity through the antisymmetric
    spectrum S^-_Φ(f). Classical systems give Δ(-ln V) = 0.

    Args:
        clock_A_filter: Dict with Ω_A(f) filter function for clock A
        clock_B_filter: Dict with Ω_B(f) filter function for clock B
        antisymmetric_psd: Dict with S^-_Φ(f) from spectral_core
        frequency_var: Frequency variable (from psd_result)
        integration_params: Optional integration bounds and method

    Returns:
        Dict with:
            'q1_witness': Δ(-ln V) value
            'integrand': Ω_A(f) Ω_B*(f) S^-_Φ(f)
            'quantum_signature': magnitude of witness
            'frequency_domain': integration details
    """
    if integration_params is None:
        integration_params = {'method': 'symbolic', 'bounds': 'infinite'}

    # Extract filter functions
    f = frequency_var
    Omega_A = clock_A_filter['filter_function']
    Omega_B = clock_B_filter['filter_function']
    S_minus_Phi = antisymmetric_psd['antisymmetric_psd']

    # Complex conjugate of Ω_B
    Omega_B_conj = conjugate(Omega_B)

    # Integrand: Ω_A(f) Ω_B*(f) S^-_Φ(f)
    integrand = Omega_A * Omega_B_conj * S_minus_Phi

    # Q1 witness: Δ(-ln V) = 2 Im ∫ df [integrand]
    if integration_params['method'] == 'symbolic':
        # Symbolic integration over infinite domain
        if integration_params['bounds'] == 'infinite':
            integral_result = integrate(integrand, f, -oo, oo)
        else:
            f_min, f_max = integration_params['bounds']
            integral_result = integrate(integrand, f, f_min, f_max)

        # Extract imaginary part and multiply by 2
        q1_witness = 2 * integral_result.imag_part()

    elif integration_params['method'] == 'series_expansion':
        # Low-frequency expansion for analytical insight
        integrand_expanded = integrand.taylor(f, 0, 6)
        # Integrate term by term (requires careful handling of convergence)
        integral_expanded = integrate(integrand_expanded, f, -oo, oo)
        q1_witness = 2 * integral_expanded.imag_part()

    else:
        raise ValueError(f"Unknown integration method: {integration_params['method']}")

    # Quantum signature magnitude
    quantum_signature = abs(q1_witness)

    return {
        'q1_witness': q1_witness,
        'integrand': integrand,
        'quantum_signature': quantum_signature,
        'frequency_domain': {
            'frequency_var': f,
            'integration_method': integration_params['method'],
            'integration_bounds': integration_params['bounds']
        },
        'filter_A': Omega_A,
        'filter_B': Omega_B_conj,
        'antisymmetric_spectrum': S_minus_Phi
    }

def compute_visibility_difference(forward_protocol, reverse_protocol):
    """
    Compute visibility difference between time-forward and time-reversed protocols.

    The quantum witness emerges from the difference:
    Δ(-ln V) = (-ln V_forward) - (-ln V_reverse)

    For classical systems: V_forward = V_reverse ⟹ Δ(-ln V) = 0
    For quantum systems: V_forward ≠ V_reverse ⟹ Δ(-ln V) ≠ 0

    Args:
        forward_protocol: Dict with visibility V_forward
        reverse_protocol: Dict with visibility V_reverse

    Returns:
        Dict with:
            'visibility_difference': Δ(-ln V)
            'forward_visibility': V_forward
            'reverse_visibility': V_reverse
            'asymmetry_parameter': |Δ(-ln V)| / <ln V>
    """
    V_forward = forward_protocol['visibility']
    V_reverse = reverse_protocol['visibility']

    # Logarithmic visibilities
    ln_V_forward = -log(V_forward)
    ln_V_reverse = -log(V_reverse)

    # Visibility difference: Δ(-ln V)
    visibility_difference = ln_V_forward - ln_V_reverse

    # Asymmetry parameter: normalized difference
    average_ln_V = (ln_V_forward + ln_V_reverse) / 2
    asymmetry_parameter = abs(visibility_difference) / average_ln_V

    return {
        'visibility_difference': visibility_difference,
        'forward_visibility': V_forward,
        'reverse_visibility': V_reverse,
        'forward_ln_visibility': ln_V_forward,
        'reverse_ln_visibility': ln_V_reverse,
        'asymmetry_parameter': asymmetry_parameter
    }

def extract_quantum_signature(q1_result, classical_comparison=None):
    """
    Extract and characterize the quantum signature from Q1 witness.

    Analyzes the magnitude, frequency dependence, and distinguishability
    of the quantum signature from classical backgrounds.

    Args:
        q1_result: Result from compute_q1_witness
        classical_comparison: Optional classical limit for comparison

    Returns:
        Dict with:
            'signature_magnitude': |Δ(-ln V)|
            'signature_phase': arg(Δ(-ln V))
            'frequency_characteristics': dominant frequency scales
            'classical_deviation': deviation from classical limit
    """
    q1_witness = q1_result['q1_witness']
    integrand = q1_result['integrand']
    f = q1_result['frequency_domain']['frequency_var']

    # Signature properties
    signature_magnitude = abs(q1_witness)
    signature_phase = arg(q1_witness)

    # Frequency characteristics: find dominant contributions
    integrand_magnitude = abs(integrand)

    # For symbolic analysis, identify key frequency scales
    # This requires examining the structure of the integrand
    frequency_characteristics = {
        'integrand_form': str(integrand),
        'dominant_scales': 'analysis_requires_specific_parameters'
    }

    # Classical deviation
    if classical_comparison is not None:
        classical_witness = classical_comparison.get('classical_q1_witness', 0)
        classical_deviation = abs(q1_witness - classical_witness)
    else:
        # Classical limit should be zero
        classical_deviation = signature_magnitude

    return {
        'signature_magnitude': signature_magnitude,
        'signature_phase': signature_phase,
        'frequency_characteristics': frequency_characteristics,
        'classical_deviation': classical_deviation,
        'quantum_vs_classical_ratio': signature_magnitude / max(abs(classical_comparison or 0), 1e-16)
    }

def verify_classical_limit(clock_setup, quantum_params):
    """
    Verify that classical systems produce zero Q1 witness.

    In the classical limit (ħ → 0), the antisymmetric spectrum S^-_Φ → 0,
    which should make the Q1 witness vanish: Δ(-ln V) → 0.

    Args:
        clock_setup: Clock configuration (filters, baselines)
        quantum_params: Quantum parameters including ħ

    Returns:
        Dict with:
            'classical_limit_witness': Q1 in ħ → 0 limit
            'limit_verified': boolean verification
            'hbar_scaling': how Q1 scales with ħ
    """
    hbar = quantum_params['hbar']

    # Build antisymmetric spectrum with explicit ħ dependence
    from .spectral_core import build_antisymmetric_spectrum, build_psd_lapse_fluctuations

    # Classical flux law PSD (no quantum corrections)
    flux_params = {
        'G': var('G'),
        'c': var('c'),
        'R': var('R'),
        'energy_density_psd': var('S_P'),
        'cutoff_frequency': var('f_c')  # Symbolic cutoff for classical limit analysis
    }
    classical_psd = build_psd_lapse_fluctuations(flux_params)

    # Antisymmetric spectrum with ħ dependence
    antisymmetric_data = build_antisymmetric_spectrum(classical_psd, quantum_params)

    # Q1 witness with explicit ħ
    clock_A = clock_setup['clock_A']
    clock_B = clock_setup['clock_B']
    q1_classical = compute_q1_witness(clock_A, clock_B, antisymmetric_data, classical_psd['frequency'])

    # Take ħ → 0 limit
    classical_limit_witness = limit(q1_classical['q1_witness'], hbar, 0)

    # Verify limit is zero
    limit_verified = bool(classical_limit_witness == 0)

    # Analyze ħ scaling
    q1_series = q1_classical['q1_witness'].taylor(hbar, 0, 3)
    hbar_scaling = q1_series

    return {
        'classical_limit_witness': classical_limit_witness,
        'limit_verified': limit_verified,
        'hbar_scaling': hbar_scaling,
        'verification_details': {
            'quantum_parameter': hbar,
            'classical_limit': 'hbar_to_zero',
            'expected_result': 0
        }
    }

def compute_time_reversal_protocol(base_protocol, reversal_params=None):
    """
    Implement time-reversal protocol for quantum witness detection.

    From Paper V: Run sequences with reversed time ordering to isolate
    the antisymmetric component of the field correlations.

    Args:
        base_protocol: Original clock sequence protocol
        reversal_params: Parameters for time reversal implementation

    Returns:
        Dict with:
            'reversed_protocol': Time-reversed sequence
            'reversal_transformation': T → -T mapping
            'phase_conjugation': complex conjugation rules
    """
    if reversal_params is None:
        reversal_params = {'method': 'phase_conjugation'}

    original_filter = base_protocol['filter_function']
    omega = base_protocol['omega']

    if reversal_params['method'] == 'phase_conjugation':
        # Time reversal: ω → -ω, then complex conjugation
        # Extract omega variable from the filter expression
        filter_variables = list(original_filter.variables())
        omega_from_filter = None
        for var in filter_variables:
            if str(var) == 'omega':
                omega_from_filter = var
                break

        if omega_from_filter is None:
            raise ValueError(f"No omega variable found in filter expression. Variables: {filter_variables}")

        # Use the actual omega variable from the filter
        reversed_filter = conjugate(original_filter.subs({omega_from_filter: -omega_from_filter}))

    elif reversal_params['method'] == 'sequence_reversal':
        # Reverse pulse sequence ordering
        # This requires specific sequence structure
        reversed_filter = original_filter  # Placeholder for sequence-specific reversal

    else:
        raise ValueError(f"Unknown reversal method: {reversal_params['method']}")

    return {
        'reversed_protocol': {
            'filter_function': reversed_filter,
            'omega': omega,
            'reversal_applied': True
        },
        'reversal_transformation': reversal_params['method'],
        'phase_conjugation': 'complex_conjugate_applied'
    }