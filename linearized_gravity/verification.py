#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Linearized Gravity Verification Script

Simple verification that Standard GR ≡ Lapse-First GR
Runs both approaches and compares results for mathematical equivalence.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

# Import our clean modules
import standard_gr
import lapse_first_gr
import equivalence_proof

def run_verification():
    """
    Complete verification of Standard GR ≡ Lapse-First GR equivalence

    Returns:
        bool: True if equivalence verified, False otherwise
    """
    try:
        # Run Standard GR computation
        print("Computing Standard GR...")
        standard_metric_data = standard_gr.compute_standard_linearization()
        standard_tt_data = standard_gr.apply_tt_gauge(standard_metric_data)
        standard_result = standard_gr.extract_physical_modes(standard_metric_data, standard_tt_data)

        # Run Lapse-First GR computation
        print("Computing Lapse-First GR...")
        lapse_variables_data = lapse_first_gr.compute_lapse_first_variables()
        lapse_constraint_data = lapse_first_gr.apply_adm_constraints(lapse_variables_data)
        lapse_first_result = lapse_first_gr.extract_tt_modes(lapse_variables_data, lapse_constraint_data)

        # Verify equivalence
        print("Verifying mathematical equivalence...")
        equivalence_result = equivalence_proof.complete_equivalence_verification(
            standard_result, lapse_first_result
        )

        # Report results
        print("\n" + "="*60)
        print("GATG LINEARIZED GRAVITY VERIFICATION RESULTS")
        print("="*60)

        dof_check = equivalence_result['dof_verification']
        print(f"DOF Equivalence: {dof_check['equivalence_verified']}")
        print(f"  Standard GR final DOF: {dof_check['standard_final_dof']}")
        print(f"  Lapse-First GR final DOF: {dof_check['lapse_first_final_dof']}")

        wave_check = equivalence_result['wave_operator_verification']
        print(f"Wave Operator Equivalence: {wave_check['wave_operator_equivalence']}")

        speed_check = equivalence_result['propagation_speed_verification']
        print(f"Propagation Speed Equivalence: {speed_check['speed_verification']}")

        overall_success = equivalence_result['mathematical_equivalence_proven']

        print("\n" + "="*60)
        if overall_success:
            print("✓ EQUIVALENCE VERIFIED: Standard GR ≡ Lapse-First GR")
            print("✓ GATG HYPOTHESIS CONFIRMED")
        else:
            print("✗ EQUIVALENCE VERIFICATION FAILED")
            print("✗ GATG HYPOTHESIS NOT CONFIRMED")
        print("="*60)

        return overall_success

    except Exception as e:
        print(f"\n✗ VERIFICATION ERROR: {e}")
        return False

if __name__ == "__main__":
    success = run_verification()
    sys.exit(0 if success else 1)