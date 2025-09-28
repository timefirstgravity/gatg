#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Gravitoelectromagnetism Verification Script

Simple verification that Standard GR ≡ Lapse-First GR for GEM formulations
Runs both approaches and compares results for mathematical equivalence.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

# Import our clean modules
import standard_gem
import lapse_first_gem
import equivalence_proof

def run_gem_verification():
    """
    Complete verification of Standard GR ≡ Lapse-First GR for GEM formulations

    Returns:
        bool: True if equivalence verified, False otherwise
    """
    try:
        # Run Standard GEM computation
        print("Computing Standard GEM...")
        standard_weak_field = standard_gem.compute_weak_field_metric()
        standard_gem_fields = standard_gem.derive_gem_fields(standard_weak_field)
        standard_field_equations = standard_gem.derive_gem_field_equations(standard_gem_fields, standard_weak_field)
        standard_forces = standard_gem.compute_gem_forces(standard_gem_fields)

        # Run Lapse-First GEM computation
        print("Computing Lapse-First GEM...")
        lapse_first_setup = lapse_first_gem.setup_lapse_first_gem_variables()
        lapse_first_gem_fields = lapse_first_gem.extract_gem_fields_from_lapse_first(lapse_first_setup)
        lapse_first_equations = lapse_first_gem.derive_lapse_first_field_equations(lapse_first_gem_fields, lapse_first_setup)
        lapse_first_forces = lapse_first_gem.compute_lapse_first_gem_forces(lapse_first_gem_fields)
        lapse_first_consistency = lapse_first_gem.verify_lapse_first_consistency(
            lapse_first_setup, lapse_first_gem_fields, lapse_first_equations
        )

        # Package results for equivalence checking
        standard_results = {
            'weak_field_setup': standard_weak_field,
            'gem_fields': standard_gem_fields,
            'field_equations': standard_field_equations,
            'forces': standard_forces
        }

        lapse_first_results = {
            'setup': lapse_first_setup,
            'gem_fields': lapse_first_gem_fields,
            'field_equations': lapse_first_equations,
            'forces': lapse_first_forces,
            'consistency': lapse_first_consistency
        }

        # Verify equivalence
        print("Verifying mathematical equivalence...")
        equivalence_result = equivalence_proof.complete_gem_equivalence_verification(
            standard_results, lapse_first_results
        )

        # Report results
        print("\n" + "="*60)
        print("GRAVITOELECTROMAGNETISM VERIFICATION RESULTS")
        print("="*60)

        field_check = equivalence_result['field_verification']
        print(f"GEM Field Equivalence: {field_check['field_definitions_compatible']}")
        print(f"  Gravitoelectric: {field_check['gravitoelectric_equivalence']}")
        print(f"  Gravitomagnetic: {field_check['gravitomagnetic_equivalence']}")

        equation_check = equivalence_result['equation_verification']
        print(f"Field Equations: {equation_check['field_equations_equivalent']}")
        print(f"  Gauss law equivalence: {equation_check['gauss_law_equivalence']}")
        print(f"  Momentum-Ampère correspondence: {equation_check['momentum_ampere_correspondence']}")

        force_check = equivalence_result['force_verification']
        print(f"Force Laws: {force_check['force_expressions_equivalent']}")
        print(f"  Lorentz-like structure: {force_check['lorentz_like_structure_verified']}")

        interpretation_check = equivalence_result['interpretation_verification']
        print(f"Physical Interpretation: {interpretation_check['physical_consistency']}")
        print(f"  Geometric complementarity: {interpretation_check['geometric_complementarity']}")

        # Lapse-first consistency
        consistency_check = lapse_first_results['consistency']
        print(f"Lapse-First Consistency: {consistency_check['lapse_first_gem_consistent']}")

        overall_success = equivalence_result['mathematical_equivalence_proven']

        print("\n" + "="*60)
        if overall_success:
            print("✓ EQUIVALENCE VERIFIED: Standard GR ≡ Lapse-First GR")
            print("✓ GRAVITOELECTROMAGNETISM FORMULATIONS CONFIRMED")
        else:
            print("✗ EQUIVALENCE VERIFICATION FAILED")
            print("✗ GRAVITOELECTROMAGNETISM FORMULATIONS NOT CONFIRMED")
        print("="*60)

        return overall_success

    except Exception as e:
        print(f"\n✗ VERIFICATION ERROR: {e}")
        return False

if __name__ == "__main__":
    success = run_gem_verification()
    sys.exit(0 if success else 1)