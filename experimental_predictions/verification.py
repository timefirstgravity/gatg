#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Experimental Predictions Verification Script

Simple verification that Standard GR ≡ Lapse-First GR for experimental predictions
Runs both approaches and compares results for mathematical equivalence.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

# Import our clean modules
import standard_predictions
import lapse_first_predictions
import equivalence_proof

def run_experimental_predictions_verification():
    """
    Complete verification of Standard GR ≡ Lapse-First GR for experimental predictions

    Returns:
        bool: True if equivalence verified, False otherwise
    """
    try:
        # Run Standard GR predictions
        print("Computing Standard GR experimental predictions...")
        standard_perihelion = standard_predictions.compute_perihelion_precession()
        standard_light_deflection = standard_predictions.compute_light_deflection()
        standard_redshift = standard_predictions.compute_gravitational_redshift()
        standard_waves = standard_predictions.compute_gravitational_waves()
        standard_time_delay = standard_predictions.compute_time_delay()

        # Run Lapse-First GR predictions
        print("Computing Lapse-First GR experimental predictions...")
        lapse_first_perihelion = lapse_first_predictions.compute_lapse_first_perihelion_precession()
        lapse_first_light_deflection = lapse_first_predictions.compute_lapse_first_light_deflection()
        lapse_first_redshift = lapse_first_predictions.compute_lapse_first_redshift()
        lapse_first_waves = lapse_first_predictions.compute_lapse_first_gravitational_waves()
        lapse_first_time_delay = lapse_first_predictions.compute_lapse_first_time_delay()
        lapse_first_consistency = lapse_first_predictions.verify_lapse_first_consistency()

        # Package results for equivalence checking
        standard_results = {
            'perihelion': standard_perihelion,
            'light_deflection': standard_light_deflection,
            'redshift': standard_redshift,
            'gravitational_waves': standard_waves,
            'time_delay': standard_time_delay
        }

        lapse_first_results = {
            'perihelion': lapse_first_perihelion,
            'light_deflection': lapse_first_light_deflection,
            'redshift': lapse_first_redshift,
            'gravitational_waves': lapse_first_waves,
            'time_delay': lapse_first_time_delay,
            'consistency': lapse_first_consistency
        }

        # Verify equivalence
        print("Verifying experimental prediction equivalence...")
        equivalence_result = equivalence_proof.complete_experimental_predictions_equivalence_verification(
            standard_results, lapse_first_results
        )

        # Report results
        print("\n" + "="*60)
        print("EXPERIMENTAL PREDICTIONS VERIFICATION RESULTS")
        print("="*60)

        # Perihelion precession
        perihelion_check = equivalence_result['perihelion_verification']
        print(f"Perihelion Precession: {perihelion_check['perihelion_equivalence_verified']}")
        print(f"  Standard: {standard_perihelion['formula']}")
        print(f"  Lapse-First: {lapse_first_perihelion['lapse_first_formula']}")

        # Light deflection
        deflection_check = equivalence_result['light_deflection_verification']
        print(f"Light Deflection: {deflection_check['light_deflection_equivalence_verified']}")
        print(f"  Standard: {standard_light_deflection['formula']}")
        print(f"  Lapse-First: {lapse_first_light_deflection['lapse_first_formula']}")

        # Gravitational redshift
        redshift_check = equivalence_result['redshift_verification']
        print(f"Gravitational Redshift: {redshift_check['redshift_equivalence_verified']}")
        print(f"  Standard: {standard_redshift['formula']}")
        print(f"  Lapse-First: {lapse_first_redshift['lapse_first_formula']}")

        # Gravitational waves
        waves_check = equivalence_result['gravitational_waves_verification']
        print(f"Gravitational Waves: {waves_check['gravitational_wave_equivalence_verified']}")
        print(f"  Strain amplitudes: {waves_check['strain_amplitudes_identical']}")
        print(f"  Wave speeds: {waves_check['wave_speeds_identical']}")
        print(f"  Energy loss: {waves_check['energy_loss_identical']}")

        # Time delay
        delay_check = equivalence_result['time_delay_verification']
        print(f"Shapiro Time Delay: {delay_check['time_delay_equivalence_verified']}")
        print(f"  Standard: {standard_time_delay['formula']}")
        print(f"  Lapse-First: {lapse_first_time_delay['lapse_first_formula']}")

        # Physical interpretations
        interpretation_check = equivalence_result['interpretation_verification']
        print(f"Physical Interpretations: {interpretation_check['interpretations_complementary']}")
        print(f"  Unified understanding: {interpretation_check['unified_understanding']}")

        # Lapse-first consistency
        consistency_check = lapse_first_results['consistency']
        print(f"Lapse-First Framework: {consistency_check['gatg_framework_valid']}")

        overall_success = equivalence_result['mathematical_equivalence_proven']

        print("\n" + "="*60)
        if overall_success:
            print("✓ EQUIVALENCE VERIFIED: Standard GR ≡ Lapse-First GR")
            print("✓ EXPERIMENTAL PREDICTIONS CONFIRMED")
            print("  • Same perihelion precession: 6πGM/(ac²(1-e²))")
            print("  • Same light deflection: 4GM/(bc²)")
            print("  • Same gravitational redshift from lapse functions")
            print("  • Same gravitational wave properties")
            print("  • Same Shapiro time delay")
        else:
            print("✗ EQUIVALENCE VERIFICATION FAILED")
            print("✗ EXPERIMENTAL PREDICTIONS NOT CONFIRMED")
        print("="*60)

        return overall_success

    except Exception as e:
        print(f"\n✗ VERIFICATION ERROR: {e}")
        return False

if __name__ == "__main__":
    success = run_experimental_predictions_verification()
    sys.exit(0 if success else 1)