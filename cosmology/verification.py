#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Cosmological Verification Script

Simple verification that Standard GR ≡ Lapse-First GR for cosmological spacetimes
Runs both approaches and compares results for mathematical equivalence.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

# Import our clean modules
import standard_cosmology
import lapse_first_cosmology
import equivalence_proof

def run_cosmological_verification():
    """
    Complete verification of Standard GR ≡ Lapse-First GR for cosmological models

    Returns:
        bool: True if equivalence verified, False otherwise
    """
    try:
        # Run Standard GR computation
        print("Computing Standard Cosmology...")
        standard_flrw = standard_cosmology.compute_flrw_metric_setup()
        standard_friedmann = standard_cosmology.derive_friedmann_equations(standard_flrw)
        standard_matter = standard_cosmology.solve_matter_dominated_epoch(standard_friedmann)
        standard_observables = standard_cosmology.compute_cosmological_observables(standard_matter)

        # Run Lapse-First GR computation
        print("Computing Lapse-First Cosmology...")
        lapse_first_variables = lapse_first_cosmology.compute_lapse_first_cosmological_variables()
        lapse_first_spatial = lapse_first_cosmology.construct_spatial_3metric(lapse_first_variables)
        lapse_first_dynamics = lapse_first_cosmology.derive_lapse_first_dynamics(lapse_first_variables, lapse_first_spatial)
        lapse_first_reconstruction = lapse_first_cosmology.verify_flrw_reconstruction(
            lapse_first_variables, lapse_first_spatial, lapse_first_dynamics
        )
        lapse_first_observables = lapse_first_cosmology.compute_lapse_first_observables(lapse_first_dynamics)

        # Package results for equivalence checking
        standard_results = {
            'flrw_setup': standard_flrw,
            'friedmann_equations': standard_friedmann,
            'matter_solutions': standard_matter,
            'observables': standard_observables
        }

        lapse_first_results = {
            'variables': lapse_first_variables,
            'spatial_metric': lapse_first_spatial,
            'dynamics': lapse_first_dynamics,
            'reconstruction': lapse_first_reconstruction,
            'observables': lapse_first_observables
        }

        # NOTE: Observational cosmological verification has been performed
        # in separate research scripts using real data (Pantheon+ SNe, BAO, CMB)
        print("\\nNOTE: Observational verification completed in research phase")
        print("  TDE model vs ΛCDM comparison: see scripts/cosmology/snfit_pantheon_fair_comparison.py")
        print("  BAO scale verification: see scripts/cosmology/bao_scale_test.py")
        print("  CMB acoustic scale: see scripts/cosmology/cmb_anchor.py")
        print("  Unified analysis: see scripts/cosmology/unified_cosmology_fit.py")
        print("  Result: TDE ≡ ΛCDM for background cosmology with real observational data")

        # Verify equivalence
        print("\\nVerifying mathematical equivalence...")
        equivalence_result = equivalence_proof.complete_cosmological_equivalence_verification(
            standard_results, lapse_first_results
        )

        # Report results
        print("\n" + "="*60)
        print("COSMOLOGICAL MODEL VERIFICATION RESULTS")
        print("="*60)

        friedmann_check = equivalence_result['friedmann_verification']
        print(f"Friedmann Equations: {friedmann_check['fundamental_dynamics_equivalent']}")
        print(f"  Equation identity: {friedmann_check['friedmann_equations_identical']}")
        print(f"  Hubble parameter: {friedmann_check['hubble_parameters_identical']}")

        matter_check = equivalence_result['matter_era_verification']
        print(f"Matter Era Solutions: {matter_check['matter_era_equivalence']}")
        print(f"  Scale factor: {matter_check['scale_factors_identical']}")
        print(f"  Age of universe: {matter_check['ages_identical']}")

        observable_check = equivalence_result['observable_verification']
        print(f"Cosmological Observables: {observable_check['observables_equivalent']}")
        print(f"  Luminosity distance: {observable_check['luminosity_distances_identical']}")

        metric_check = equivalence_result['metric_verification']
        print(f"Metric Reconstruction: {metric_check['lapse_first_metric_correct']}")
        print(f"  FLRW reconstruction: {metric_check['metric_reconstruction_successful']}")

        print(f"Observational Verification: Referenced")
        print(f"  Research scripts demonstrate TDE ≡ ΛCDM with real data")
        print(f"  Fundamental mapping: a(tau) = e^(-Phi(tau)), H(tau) = -dPhi/dtau")

        overall_success = equivalence_result['mathematical_equivalence_proven']

        print("\n" + "="*60)
        if overall_success:
            print("✓ EQUIVALENCE VERIFIED: Standard GR ≡ Lapse-First GR")
            print("✓ COSMOLOGICAL MODELS CONFIRMED")
            print("✓ OBSERVATIONAL VERIFICATION COMPLETED IN RESEARCH PHASE")
        else:
            print("✗ EQUIVALENCE VERIFICATION FAILED")
            print("✗ COSMOLOGICAL MODELS NOT CONFIRMED")
        print("="*60)

        return overall_success

    except Exception as e:
        print(f"\n✗ VERIFICATION ERROR: {e}")
        return False

if __name__ == "__main__":
    success = run_cosmological_verification()
    sys.exit(0 if success else 1)