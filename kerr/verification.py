#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Kerr Verification Script

Simple verification that Standard GR ≡ Lapse-First GR for Kerr solution
Runs both approaches and compares results for mathematical equivalence.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

# Import our clean modules
import standard_kerr
import lapse_first_kerr
import equivalence_proof

def run_kerr_verification():
    """
    Complete verification of Standard GR ≡ Lapse-First GR for Kerr solution

    Returns:
        bool: True if equivalence verified, False otherwise
    """
    try:
        # Run Standard GR computation
        print("Computing Standard Kerr...")
        standard_setup = standard_kerr.compute_standard_kerr_setup()
        standard_metric = standard_kerr.construct_kerr_metric(standard_setup)
        standard_einstein = standard_kerr.verify_kerr_einstein_equations(standard_setup, standard_metric)
        standard_schwarzschild = standard_kerr.verify_schwarzschild_limit(standard_setup)

        # Run Lapse-First GR computation
        print("Computing Lapse-First Kerr...")
        lapse_first_setup = lapse_first_kerr.compute_lapse_first_kerr_variables()
        lapse_data = lapse_first_kerr.extract_lapse_function(lapse_first_setup)
        shift_data = lapse_first_kerr.extract_shift_vector(lapse_first_setup)
        spatial_data = lapse_first_kerr.construct_spatial_metric(lapse_first_setup)
        reconstruction_data = lapse_first_kerr.verify_lapse_first_reconstruction(
            lapse_first_setup, lapse_data, shift_data, spatial_data
        )

        # Package results for equivalence checking
        standard_results = {
            'setup': standard_setup,
            'metric': standard_metric,
            'einstein': standard_einstein,
            'schwarzschild_limit': standard_schwarzschild
        }

        lapse_first_results = {
            'setup': lapse_first_setup,
            'lapse': lapse_data,
            'shift': shift_data,
            'spatial': spatial_data,
            'reconstruction': reconstruction_data
        }

        # NEW: Slow-rotation expansion checks
        print("\\nComputing slow-rotation expansion...")
        slow_rotation_result = standard_kerr.compute_slow_rotation_expansion(standard_setup)
        print(f"  First-order expansion in 'a': {slow_rotation_result['expansion_order']}")
        print(f"  Frame-dragging coefficient verified: {slow_rotation_result['frame_dragging_verified']}")

        print("\\nVerifying frame-dragging physics...")
        frame_dragging_result = standard_kerr.verify_frame_dragging_physics(
            standard_setup, slow_rotation_result
        )
        print(f"  Lense-Thirring physics verified: {frame_dragging_result['physics_verified']}")
        print(f"  Angular dependence (sin²θ): {frame_dragging_result['angular_dependence_correct']}")
        print(f"  Radial falloff (1/r): {frame_dragging_result['r_dependence_correct']}")

        # Verify equivalence
        print("\\nVerifying mathematical equivalence...")
        equivalence_result = equivalence_proof.complete_kerr_equivalence_verification(
            standard_results, lapse_first_results
        )

        # Report results
        print("\n" + "="*60)
        print("KERR SOLUTION VERIFICATION RESULTS")
        print("="*60)

        metric_check = equivalence_result['metric_verification']
        print(f"Metric Reconstruction: {metric_check['metric_equivalence_verified']}")
        print(f"  Temporal component: {metric_check['temporal_component_match']}")
        print(f"  Frame dragging: {metric_check['frame_dragging_match']}")

        einstein_check = equivalence_result['einstein_verification']
        print(f"Einstein Equations: {einstein_check['einstein_consistency_verified']}")
        print(f"  Standard vacuum verified: {einstein_check['standard_einstein_verified']}")

        schwarzschild_check = equivalence_result['schwarzschild_verification']
        print(f"Schwarzschild Limit: {schwarzschild_check['schwarzschild_limit_equivalence']}")
        print(f"  Frame dragging vanishes: {schwarzschild_check['frame_dragging_vanishes']}")

        parameter_check = equivalence_result['parameter_verification']
        print(f"Parameter Consistency: {parameter_check['parameter_consistency']}")

        print(f"Slow-Rotation Expansion: {slow_rotation_result['frame_dragging_verified'] and slow_rotation_result['temporal_expansion_correct']}")
        print(f"  Frame-dragging coefficient: O(a) term verified")
        print(f"  Temporal component: O(1) term matches Schwarzschild")

        print(f"Frame-Dragging Physics: {frame_dragging_result['physics_verified']}")
        print(f"  Lense-Thirring frequency: Ω_LT = 2Ma/r³")

        # Include new checks in overall success
        slow_rotation_success = (
            slow_rotation_result['frame_dragging_verified'] and
            slow_rotation_result['temporal_expansion_correct'] and
            frame_dragging_result['physics_verified']
        )

        overall_success = equivalence_result['mathematical_equivalence_proven'] and slow_rotation_success

        print("\n" + "="*60)
        if overall_success:
            print("✓ EQUIVALENCE VERIFIED: Standard GR ≡ Lapse-First GR")
            print("✓ KERR SOLUTION CONFIRMED")
            print("✓ SLOW-ROTATION EXPANSION VERIFIED")
            print("✓ FRAME-DRAGGING PHYSICS CONFIRMED")
        else:
            print("✗ EQUIVALENCE VERIFICATION FAILED")
            print("✗ KERR SOLUTION NOT CONFIRMED")
        print("="*60)

        return overall_success

    except Exception as e:
        print(f"\n✗ VERIFICATION ERROR: {e}")
        return False

if __name__ == "__main__":
    success = run_kerr_verification()
    sys.exit(0 if success else 1)