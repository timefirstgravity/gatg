#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Schwarzschild Verification Script

Simple verification that Standard GR ≡ Lapse-First GR for Schwarzschild solution
Runs both approaches and compares results for mathematical equivalence.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

# Import our clean modules
import standard_schwarzschild
import lapse_first_schwarzschild
import equivalence_proof

def run_schwarzschild_verification():
    """
    Complete verification of Standard GR ≡ Lapse-First GR for Schwarzschild solution

    Returns:
        bool: True if equivalence verified, False otherwise
    """
    try:
        # Run Standard GR computation
        print("Computing Standard Schwarzschild...")
        standard_setup = standard_schwarzschild.compute_standard_schwarzschild_setup()
        standard_ode = standard_schwarzschild.derive_schwarzschild_ode(standard_setup)
        standard_solution = standard_schwarzschild.solve_schwarzschild_ode(standard_ode)
        standard_result = standard_schwarzschild.verify_schwarzschild_solution(standard_setup, standard_solution)

        # Run Lapse-First GR computation
        print("Computing Lapse-First Schwarzschild...")
        lapse_variables = lapse_first_schwarzschild.compute_lapse_first_variables()
        lapse_spatial = lapse_first_schwarzschild.construct_spatial_metric(lapse_variables)
        lapse_metric = lapse_first_schwarzschild.construct_lapse_first_metric(lapse_variables, lapse_spatial)
        lapse_first_result = lapse_first_schwarzschild.verify_lapse_equivalence(lapse_variables, lapse_metric)

        # Combine results for equivalence checking
        standard_combined = {**standard_result, **standard_solution}
        lapse_first_combined = {**lapse_first_result, **lapse_variables, **lapse_metric}

        # Verify equivalence
        print("Verifying mathematical equivalence...")
        equivalence_result = equivalence_proof.complete_schwarzschild_equivalence_verification(
            standard_combined, lapse_first_combined
        )

        # NEW: Compute Kretschmann scalar
        print("\nComputing Kretschmann scalar K = R_μνρσ R^μνρσ...")
        kretschmann_result = standard_schwarzschild.compute_kretschmann_scalar(
            standard_setup, standard_result
        )
        print(f"  Kretschmann scalar: K = {kretschmann_result['kretschmann_scalar']}")
        print(f"  Expected value:     K = {kretschmann_result['expected_value']}")
        print(f"  Match verified: {kretschmann_result['match_verified']}")

        # NEW: Verify Birkhoff's theorem
        print("\nVerifying Birkhoff's theorem (uniqueness)...")
        birkhoff_result = standard_schwarzschild.verify_birkhoff_theorem(
            standard_ode, standard_solution
        )
        print(f"  Birkhoff theorem verified: {birkhoff_result['birkhoff_verified']}")
        print(f"  General solution: A(r) = {birkhoff_result['general_solution']}")
        print(f"  Unique physical solution: A(r) = {birkhoff_result['unique_physical_solution']}")
        print(f"  Uniqueness condition: {birkhoff_result['uniqueness_condition']}")

        # Report results
        print("\n" + "="*60)
        print("SCHWARZSCHILD SOLUTION VERIFICATION RESULTS")
        print("="*60)

        metric_check = equivalence_result['metric_verification']
        print(f"Metric Equivalence: {metric_check['metrics_identical']}")

        ode_check = equivalence_result['ode_verification']
        print(f"ODE Consistency: {ode_check['ode_consistency_verified']}")
        print(f"  Standard A(r): {ode_check['standard_A']}")
        print(f"  Lapse-First A(r): {ode_check['lapse_first_A']}")

        vacuum_check = equivalence_result['vacuum_verification']
        print(f"Vacuum Conditions: {vacuum_check['vacuum_equivalence_proven']}")

        # Include new verifications in overall success
        overall_success = (
            equivalence_result['mathematical_equivalence_proven'] and
            kretschmann_result['match_verified'] and
            birkhoff_result['birkhoff_verified']
        )

        print("Kretschmann Scalar: {0}".format("✓ Matches K = 12r_s²/r⁶" if kretschmann_result['match_verified'] else "✗ Does not match"))
        print("Birkhoff Theorem: {0}".format("✓ Uniqueness verified" if birkhoff_result['birkhoff_verified'] else "✗ Not verified"))

        print("\n" + "="*60)
        if overall_success:
            print("✓ EQUIVALENCE VERIFIED: Standard GR ≡ Lapse-First GR")
            print("✓ SCHWARZSCHILD SOLUTION CONFIRMED")
            print("✓ KRETSCHMANN SCALAR MATCHES")
            print("✓ BIRKHOFF UNIQUENESS THEOREM VERIFIED")
        else:
            print("✗ EQUIVALENCE VERIFICATION FAILED")
            print("✗ SCHWARZSCHILD SOLUTION NOT CONFIRMED")
        print("="*60)

        return overall_success

    except Exception as e:
        print(f"\n✗ VERIFICATION ERROR: {e}")
        return False

if __name__ == "__main__":
    success = run_schwarzschild_verification()
    sys.exit(0 if success else 1)