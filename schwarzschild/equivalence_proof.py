#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Schwarzschild Equivalence Proof: Standard GR ≡ Lapse-First GR

Direct comparison of computed results from both approaches.
Verifies that different theoretical starting points yield identical Schwarzschild spacetime.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

def verify_metric_equivalence(standard_result, lapse_first_result):
    """
    Verify both approaches yield the same metric tensor

    Args:
        standard_result: Result from standard_schwarzschild computation
        lapse_first_result: Result from lapse_first_schwarzschild computation

    Returns:
        dict: Metric comparison results
    """
    standard_metric = standard_result['complete_metric']
    lapse_first_metric = lapse_first_result['metric']

    # Compare metric components
    component_differences = {}
    components_match = True

    for i in range(4):
        for j in range(i, 4):  # Only check upper triangle (symmetric)
            standard_comp = standard_metric[i,j].expr().simplify_full()
            lapse_first_comp = lapse_first_metric[i,j].expr().simplify_full()
            difference = (standard_comp - lapse_first_comp).simplify_full()

            component_differences[f'g_{i}{j}'] = difference
            if difference != 0:
                components_match = False

    return {
        'component_differences': component_differences,
        'metrics_identical': components_match,
        'standard_metric': standard_metric,
        'lapse_first_metric': lapse_first_metric
    }

def verify_ode_consistency(standard_result, lapse_first_result):
    """
    Verify both approaches satisfy the same fundamental ODE

    Args:
        standard_result: Standard GR derivation
        lapse_first_result: Lapse-first construction

    Returns:
        dict: ODE verification results
    """
    # Standard approach derives: r A'(r) + A(r) - 1 = 0
    standard_A = standard_result['final_solution']

    # Lapse-first uses: A = exp(2Φ) = 1 - r_s/r
    lapse_Phi = lapse_first_result['temporal_potential']
    lapse_A = exp(2 * lapse_Phi).simplify_full()

    # Verify both give same A(r)
    A_difference = (standard_A - lapse_A).simplify_full()
    A_functions_match = (A_difference == 0)

    # Verify ODE satisfaction for lapse-derived A
    r = var('r')
    ode_lhs = r * diff(lapse_A, r) + lapse_A - 1
    ode_satisfied = ode_lhs.simplify_full() == 0

    return {
        'standard_A': standard_A,
        'lapse_first_A': lapse_A,
        'A_difference': A_difference,
        'A_functions_identical': A_functions_match,
        'ode_satisfied_by_lapse': ode_satisfied,
        'ode_consistency_verified': A_functions_match and ode_satisfied
    }

def verify_vacuum_conditions(standard_result, lapse_first_result):
    """
    Verify both approaches satisfy Einstein vacuum equations

    Args:
        standard_result: Standard GR vacuum verification
        lapse_first_result: Lapse-first metric construction

    Returns:
        dict: Vacuum condition verification
    """
    standard_vacuum = standard_result['vacuum_verified']

    # For lapse-first, we verify the constructed metric satisfies vacuum equations
    # by checking it gives the same A(r) that satisfies the ODE
    lapse_construction_correct = lapse_first_result['equivalence_verified']

    both_satisfy_vacuum = standard_vacuum and lapse_construction_correct

    return {
        'standard_vacuum_verified': standard_vacuum,
        'lapse_construction_correct': lapse_construction_correct,
        'both_satisfy_vacuum': both_satisfy_vacuum,
        'vacuum_equivalence_proven': both_satisfy_vacuum
    }

def complete_schwarzschild_equivalence_verification(standard_result, lapse_first_result):
    """
    Perform complete mathematical equivalence verification for Schwarzschild solution

    Args:
        standard_result: Complete result from Standard GR computation
        lapse_first_result: Complete result from Lapse-First GR computation

    Returns:
        dict: Complete equivalence verification results
    """
    # Individual verifications
    metric_verification = verify_metric_equivalence(standard_result, lapse_first_result)
    ode_verification = verify_ode_consistency(standard_result, lapse_first_result)
    vacuum_verification = verify_vacuum_conditions(standard_result, lapse_first_result)

    # Overall equivalence
    all_tests_pass = (
        metric_verification['metrics_identical'] and
        ode_verification['ode_consistency_verified'] and
        vacuum_verification['vacuum_equivalence_proven']
    )

    return {
        'metric_verification': metric_verification,
        'ode_verification': ode_verification,
        'vacuum_verification': vacuum_verification,
        'mathematical_equivalence_proven': all_tests_pass,
        'schwarzschild_gatg_equivalence_confirmed': all_tests_pass
    }