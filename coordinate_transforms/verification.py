#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Coordinate Transform Verification Script

Simple verification that Standard GR ≡ Lapse-First GR for coordinate systems
Runs both approaches and compares results for mathematical equivalence.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

# Import our clean modules
import standard_coordinates
import lapse_first_coordinates
import equivalence_proof

def run_coordinate_transform_verification():
    """
    Complete verification of Standard GR ≡ Lapse-First GR for coordinate systems

    Returns:
        bool: True if equivalence verified, False otherwise
    """
    try:
        # Run Standard coordinate computations
        print("Computing Standard coordinate systems...")
        standard_schwarzschild = standard_coordinates.setup_schwarzschild_coordinates()
        standard_ef = standard_coordinates.setup_eddington_finkelstein_coordinates()
        standard_kruskal = standard_coordinates.setup_kruskal_szekeres_coordinates()
        standard_transform = standard_coordinates.compute_coordinate_transformation_schwarzschild_to_ef()
        standard_consistency = standard_coordinates.verify_coordinate_transformation_consistency(
            standard_schwarzschild, standard_ef, standard_transform
        )

        # Run Lapse-First coordinate computations
        print("Computing Lapse-First coordinate systems...")
        lapse_first_schwarzschild = lapse_first_coordinates.setup_lapse_first_schwarzschild_form()
        lapse_first_ef = lapse_first_coordinates.setup_lapse_first_eddington_finkelstein_form()
        lapse_first_kruskal = lapse_first_coordinates.setup_lapse_first_kruskal_form()
        lapse_first_transform = lapse_first_coordinates.compute_lapse_first_coordinate_transformations()

        # Verify reconstructions
        schwarzschild_reconstruction = lapse_first_coordinates.verify_lapse_first_reconstruction(
            lapse_first_schwarzschild, 'Lapse-First Schwarzschild'
        )
        ef_reconstruction = lapse_first_coordinates.verify_lapse_first_reconstruction(
            lapse_first_ef, 'Lapse-First Eddington-Finkelstein'
        )

        # Package results for equivalence checking
        standard_results = {
            'schwarzschild': standard_schwarzschild,
            'ef': standard_ef,
            'kruskal': standard_kruskal,
            'transformation': standard_transform,
            'consistency': standard_consistency
        }

        lapse_first_results = {
            'schwarzschild': lapse_first_schwarzschild,
            'ef': lapse_first_ef,
            'kruskal': lapse_first_kruskal,
            'transformation': lapse_first_transform,
            'schwarzschild_reconstruction': schwarzschild_reconstruction,
            'ef_reconstruction': ef_reconstruction
        }

        # Verify equivalence
        print("Verifying mathematical equivalence...")
        equivalence_result = equivalence_proof.complete_coordinate_equivalence_verification(
            standard_results, lapse_first_results
        )

        # Report results
        print("\n" + "="*60)
        print("COORDINATE TRANSFORM VERIFICATION RESULTS")
        print("="*60)

        schwarzschild_check = equivalence_result['schwarzschild_verification']
        print(f"Schwarzschild Coordinates: {schwarzschild_check['schwarzschild_coordinates_equivalent']}")
        print(f"  Metric components: g_tt={schwarzschild_check['g_tt_equivalent']}, g_rr={schwarzschild_check['g_rr_equivalent']}")
        print(f"  Static spacetime: {schwarzschild_check['shift_vector_zero']}")

        ef_check = equivalence_result['eddington_finkelstein_verification']
        print(f"Eddington-Finkelstein: {ef_check['eddington_finkelstein_equivalent']}")
        print(f"  Mass consistency: {ef_check['mass_consistency']}")
        print(f"  Non-static frame: {ef_check['radial_shift_present']}")
        print(f"  Mixed components: {ef_check['mixed_component_present']}")

        transform_check = equivalence_result['transformation_verification']
        print(f"Coordinate Transformations: {transform_check['transformation_approaches_consistent']}")
        print(f"  Transform equivalence: {transform_check['coordinate_transformations_equivalent']}")

        physical_check = equivalence_result['physical_verification']
        print(f"Physical Spacetime: {physical_check['physical_spacetime_equivalent']}")
        print(f"  Mass parameters: {physical_check['mass_parameters_consistent']}")

        # Reconstruction checks
        print(f"Lapse-First Reconstructions:")
        print(f"  Schwarzschild: {schwarzschild_reconstruction['metric_reconstruction_valid']}")
        print(f"  Eddington-Finkelstein: {ef_reconstruction['metric_reconstruction_valid']}")

        overall_success = equivalence_result['mathematical_equivalence_proven']

        print("\n" + "="*60)
        if overall_success:
            print("✓ EQUIVALENCE VERIFIED: Standard GR ≡ Lapse-First GR")
            print("✓ COORDINATE SYSTEMS CONFIRMED")
        else:
            print("✗ EQUIVALENCE VERIFICATION FAILED")
            print("✗ COORDINATE SYSTEMS NOT CONFIRMED")
        print("="*60)

        return overall_success

    except Exception as e:
        print(f"\n✗ VERIFICATION ERROR: {e}")
        return False

if __name__ == "__main__":
    success = run_coordinate_transform_verification()
    sys.exit(0 if success else 1)