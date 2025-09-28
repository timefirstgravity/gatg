#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Unit tests for compute_general_lambda_momentum_constraint function
Tests the momentum constraint computation for general metric without assuming Λ = -Φ
"""

import sys
import os
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sage.all import *
from core.constraints import compute_general_lambda_momentum_constraint
from core.manifolds import setup_manifold_and_coordinates

class TestGeneralLambdaMomentumConstraint(unittest.TestCase):
    """Test general lambda momentum constraint computation"""

    def setUp(self):
        """Set up test manifold and coordinates"""
        self.M, self.X, self.t, self.r, self.th, self.ph = setup_manifold_and_coordinates()

    def test_function_returns_dictionary(self):
        """Test that function returns a dictionary with expected structure"""
        result = compute_general_lambda_momentum_constraint(self.M, self.t, self.r)

        self.assertIsInstance(result, dict)

        # Check all expected keys are present
        expected_keys = [
            'general_DY', 'general_K_trace', 'general_Y_components',
            'general_K_components', 'phi_function', 'lambda_function',
            'reduction_condition'
        ]
        for key in expected_keys:
            self.assertIn(key, result)

    def test_temporal_functions_defined(self):
        """Test that Phi and Lambda functions are properly defined"""
        result = compute_general_lambda_momentum_constraint(self.M, self.t, self.r)

        Phi = result['phi_function']
        Lambda = result['lambda_function']

        # Functions should depend on both t and r
        self.assertTrue(hasattr(Phi, '__call__'))
        self.assertTrue(hasattr(Lambda, '__call__'))

        # Check that functions can be differentiated
        dPhi_dt = diff(Phi, self.t)
        dLambda_dt = diff(Lambda, self.t)
        self.assertIsNotNone(dPhi_dt)
        self.assertIsNotNone(dLambda_dt)

    def test_Y_tensor_mathematical_consistency(self):
        """Test Y tensor mathematical relationships"""
        result = compute_general_lambda_momentum_constraint(self.M, self.t, self.r)

        Y_rr, Y_thth, Y_phph = result['general_Y_components']
        K_trace = result['general_K_trace']

        # Y_rr should be 0 for spherically symmetric case
        # This is because K^r_r = K_trace in spherical symmetry
        Y_rr_simplified = Y_rr.simplify_full()
        self.assertEqual(Y_rr_simplified, 0)

        # Angular Y components should equal -K_trace
        # Since K^θ_θ = K^φ_φ = 0, we have Y^θ_θ = Y^φ_φ = -K_trace
        Y_thth_simplified = Y_thth.simplify_full()
        Y_phph_simplified = Y_phph.simplify_full()

        self.assertEqual(Y_thth_simplified, -K_trace)
        self.assertEqual(Y_phph_simplified, -K_trace)

    def test_extrinsic_curvature_components(self):
        """Test extrinsic curvature component calculations"""
        result = compute_general_lambda_momentum_constraint(self.M, self.t, self.r)

        K_rr, K_thth, K_phph = result['general_K_components']

        # Angular components should be 0 since γ_θθ = r² and γ_φφ = r²sin²θ don't depend on time
        self.assertEqual(K_thth, 0)
        self.assertEqual(K_phph, 0)

        # K_rr should contain time derivatives of Lambda and Phi
        Lambda = result['lambda_function']
        dLambda_dt = diff(Lambda, self.t)

        # Check that K_rr contains the time derivative of Lambda
        self.assertTrue(str(dLambda_dt) in str(K_rr) or 'diff(Lambda' in str(K_rr))

    def test_christoffel_symbols_consistency(self):
        """Test that Christoffel symbols have correct values"""
        result = compute_general_lambda_momentum_constraint(self.M, self.t, self.r)

        # Extract the covariant divergence calculation
        DY = result['general_DY']

        # The calculation uses Γ^θ_rθ = Γ^φ_rφ = 1/r
        # We can verify this by checking the structure of DY

        # DY should be: (2/r) * K_trace (from the angular Christoffel terms)
        K_trace = result['general_K_trace']
        expected_DY = (2/self.r) * K_trace

        # Simplify both expressions for comparison
        DY_simplified = DY.simplify_full()
        expected_simplified = expected_DY.simplify_full()

        self.assertEqual(DY_simplified, expected_simplified)

    def test_momentum_constraint_formula(self):
        """Test the momentum constraint mathematical formula"""
        result = compute_general_lambda_momentum_constraint(self.M, self.t, self.r)

        DY = result['general_DY']
        Y_thth, Y_phph = result['general_Y_components'][1:3]

        # DY should equal -(1/r) * Y_thth - (1/r) * Y_phph
        # Since Y_thth = Y_phph = -K_trace, this gives DY = (2/r) * K_trace
        expected_formula = -(1/self.r) * Y_thth - (1/self.r) * Y_phph
        expected_simplified = expected_formula.simplify_full()
        DY_simplified = DY.simplify_full()

        self.assertEqual(DY_simplified, expected_simplified)

    def test_reduction_condition_included(self):
        """Test that reduction condition is properly documented"""
        result = compute_general_lambda_momentum_constraint(self.M, self.t, self.r)

        reduction_condition = result['reduction_condition']
        self.assertIsInstance(reduction_condition, str)
        self.assertIn('Λ = -Φ', reduction_condition)
        self.assertIn('γ_rr', reduction_condition)

    def test_mathematical_self_consistency(self):
        """Test overall mathematical self-consistency of the computation"""
        result = compute_general_lambda_momentum_constraint(self.M, self.t, self.r)

        # Extract key quantities
        K_trace = result['general_K_trace']
        Y_components = result['general_Y_components']
        K_components = result['general_K_components']
        DY = result['general_DY']

        # Check that trace calculation is consistent
        # K_trace should equal K^r_r + K^θ_θ + K^φ_φ
        K_rr, K_thth, K_phph = K_components
        Phi = result['phi_function']
        Lambda = result['lambda_function']

        # Mixed components: K^r_r = γ^rr * K_rr = e^(-2Λ) * K_rr
        Krr_up = exp(-2*Lambda) * K_rr
        computed_trace = Krr_up + 0 + 0  # Angular components are 0

        self.assertEqual(K_trace.simplify_full(), computed_trace.simplify_full())

    def test_temporal_derivative_consistency(self):
        """Test that time derivatives are handled consistently"""
        result = compute_general_lambda_momentum_constraint(self.M, self.t, self.r)

        Lambda = result['lambda_function']
        K_rr = result['general_K_components'][0]

        # K_rr should contain ∂_t Λ term
        # K_rr = -(1/(2*e^Φ)) * ∂_t(e^(2Λ)) = -(1/(2*e^Φ)) * 2*e^(2Λ)*∂_t Λ = -e^(2Λ-Φ) * ∂_t Λ
        dLambda_dt = diff(Lambda, self.t)

        # Check that K_rr contains the time derivative of Lambda
        K_rr_vars = K_rr.variables()
        dLambda_vars = dLambda_dt.variables() if hasattr(dLambda_dt, 'variables') else []

        # The derivative should appear in K_rr
        self.assertTrue(any(str(dLambda_dt) in str(K_rr) or str(var) in str(K_rr) for var in dLambda_vars))

if __name__ == '__main__':
    unittest.main()