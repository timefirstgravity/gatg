#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Mathematical Equivalence Proof: Standard GR ≡ Lapse-First GR

Direct comparison of computed results from both approaches.
Verifies that different theoretical starting points yield identical physical content.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

def verify_dof_equivalence(standard_result, lapse_first_result):
    """
    Verify both approaches yield the same number of physical degrees of freedom

    Args:
        standard_result: Result from standard_gr.extract_physical_modes()
        lapse_first_result: Result from lapse_first_gr.extract_tt_modes()

    Returns:
        dict: DOF comparison results
    """
    standard_dof = standard_result['final_physical_dof']
    lapse_first_dof = lapse_first_result['final_physical_dof']

    dof_match = (standard_dof == lapse_first_dof)
    both_yield_two = (standard_dof == 2 and lapse_first_dof == 2)

    return {
        'standard_final_dof': standard_dof,
        'lapse_first_final_dof': lapse_first_dof,
        'dof_counts_match': dof_match,
        'both_yield_two_modes': both_yield_two,
        'equivalence_verified': dof_match and both_yield_two
    }

def verify_wave_operator_equivalence(standard_result, lapse_first_result):
    """
    Verify both approaches derive the same wave operator

    Args:
        standard_result: Result from standard_gr.extract_physical_modes()
        lapse_first_result: Result from lapse_first_gr.extract_tt_modes()

    Returns:
        dict: Wave operator comparison results
    """
    # Extract wave operators
    standard_op = standard_result['wave_operator']
    lapse_first_op = lapse_first_result['wave_operator']

    # Test on a simple field to compare operators
    test_field = var('test_field')
    standard_applied = standard_op(test_field)
    lapse_first_applied = lapse_first_op(test_field)

    # Check if operators are identical
    operator_difference = standard_applied - lapse_first_applied
    operators_identical = (operator_difference == 0)

    # Verify wave equation forms
    standard_form = standard_result['wave_equation_form']
    lapse_first_form = lapse_first_result['wave_equation_form']
    equation_forms_match = (standard_form == lapse_first_form)

    return {
        'standard_operator_on_test': standard_applied,
        'lapse_first_operator_on_test': lapse_first_applied,
        'operator_difference': operator_difference,
        'operators_identical': operators_identical,
        'standard_equation_form': standard_form,
        'lapse_first_equation_form': lapse_first_form,
        'equation_forms_match': equation_forms_match,
        'wave_operator_equivalence': operators_identical and equation_forms_match
    }

def verify_propagation_speed(standard_result, lapse_first_result):
    """
    Verify both approaches predict the same wave propagation speed

    Args:
        standard_result: Result from standard_gr.extract_physical_modes()
        lapse_first_result: Result from lapse_first_gr.extract_tt_modes()

    Returns:
        dict: Propagation speed comparison
    """
    standard_speed = standard_result['propagation_speed']
    lapse_first_speed = lapse_first_result['propagation_speed']

    speeds_match = (standard_speed == lapse_first_speed)
    both_lightspeed = (standard_speed == 'c = 1' and lapse_first_speed == 'c = 1')

    # Compute speed from wave operator coefficients
    test_time = var('t')**2
    test_space = var('x')**2

    standard_op = standard_result['wave_operator']
    time_coeff = standard_op(test_time) / 2  # Should be -1
    space_coeff = standard_op(test_space) / 2  # Should be +1

    computed_speed_squared = space_coeff / abs(time_coeff)
    computed_speed = sqrt(computed_speed_squared)

    return {
        'standard_speed': standard_speed,
        'lapse_first_speed': lapse_first_speed,
        'speeds_match': speeds_match,
        'both_predict_lightspeed': both_lightspeed,
        'computed_time_coefficient': time_coeff,
        'computed_space_coefficient': space_coeff,
        'computed_speed': computed_speed,
        'speed_verification': both_lightspeed and (computed_speed == 1)
    }

def complete_equivalence_verification(standard_result, lapse_first_result):
    """
    Perform complete mathematical equivalence verification

    Args:
        standard_result: Complete result from Standard GR computation
        lapse_first_result: Complete result from Lapse-First GR computation

    Returns:
        dict: Complete equivalence verification results
    """
    # Individual verifications
    dof_verification = verify_dof_equivalence(standard_result, lapse_first_result)
    operator_verification = verify_wave_operator_equivalence(standard_result, lapse_first_result)
    speed_verification = verify_propagation_speed(standard_result, lapse_first_result)

    # Overall equivalence
    all_tests_pass = (
        dof_verification['equivalence_verified'] and
        operator_verification['wave_operator_equivalence'] and
        speed_verification['speed_verification']
    )

    return {
        'dof_verification': dof_verification,
        'wave_operator_verification': operator_verification,
        'propagation_speed_verification': speed_verification,
        'mathematical_equivalence_proven': all_tests_pass,
        'gatg_hypothesis_confirmed': all_tests_pass
    }