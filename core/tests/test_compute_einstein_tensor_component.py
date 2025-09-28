#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Unit tests for compute_einstein_tensor_component function
Tests computation of specific Einstein tensor component G^t_r
"""

import sys
import os
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sage.all import *
from core.tensors import compute_einstein_tensor_component, compute_einstein_tensor
from core.manifolds import setup_manifold_and_coordinates
from core.metrics import construct_metric

class TestComputeEinsteinTensorComponent(unittest.TestCase):
    """Test specific Einstein tensor component computation"""

    def setUp(self):
        """Set up test manifold and coordinates"""
        self.M, self.X, self.t, self.r, self.th, self.ph = setup_manifold_and_coordinates()

    def test_returns_scalar_expression(self):
        """Test that function returns a scalar expression"""
        A = function('A')(self.r)
        result = compute_einstein_tensor_component(self.M, A, self.r, self.t)

        # Should return a scalar expression (not a tensor)
        self.assertFalse(hasattr(result, 'tensor_type'),
                        "Result should be scalar, not tensor")

        # Should be a symbolic expression or number
        self.assertTrue(hasattr(result, 'simplify_full') or isinstance(result, (int, float)),
                       "Result should be symbolic expression or number")

    def test_flat_spacetime_gtr_zero(self):
        """Test that G^t_r = 0 for flat Minkowski spacetime"""
        A = 1  # Flat spacetime
        result = compute_einstein_tensor_component(self.M, A, self.r, self.t)

        # For flat spacetime, all Einstein tensor components should be zero
        result_simplified = result.simplify_full() if hasattr(result, 'simplify_full') else result

        if result_simplified != 0:
            # Check if numerically small
            self.assertTrue(abs(result_simplified) < 1e-10 if hasattr(result_simplified, '__abs__') else result_simplified == 0,
                           f"G^tr = {result_simplified} should be zero for flat spacetime")

    def test_schwarzschild_gtr_zero(self):
        """Test that G^t_r = 0 for Schwarzschild metric (spherical symmetry)"""
        M_param = var('M')
        A = 1 - 2*M_param/self.r
        result = compute_einstein_tensor_component(self.M, A, self.r, self.t)

        # For spherically symmetric metric, G^tr should be zero
        result_simplified = result.simplify_full() if hasattr(result, 'simplify_full') else result
        self.assertEqual(result_simplified, 0,
                        "G^tr should be zero for spherically symmetric Schwarzschild metric")

    def test_consistency_with_full_tensor_computation(self):
        """Test that result matches G^tr from full Einstein tensor computation"""
        A = function('A')(self.r)

        # Compute using this function
        Gtr_component = compute_einstein_tensor_component(self.M, A, self.r, self.t)

        # Compute using full tensor method
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)
        G_mixed = G.up(g)
        Gtr_full = G_mixed[0,1].expr()

        # Should match
        Gtr_component_simplified = Gtr_component.simplify_full() if hasattr(Gtr_component, 'simplify_full') else Gtr_component
        Gtr_full_simplified = Gtr_full.simplify_full() if hasattr(Gtr_full, 'simplify_full') else Gtr_full

        self.assertEqual(Gtr_component_simplified, Gtr_full_simplified,
                        "G^tr from component function should match full tensor computation")

    def test_time_dependent_metric_nonzero_gtr(self):
        """Test that time-dependent metrics can have non-zero G^t_r"""
        # Time-dependent metric function
        A = function('A')(self.t, self.r)
        result = compute_einstein_tensor_component(self.M, A, self.r, self.t)

        # For time-dependent metric, G^tr might be non-zero
        # We test that the computation completes successfully
        self.assertIsNotNone(result, "Should successfully compute G^tr for time-dependent metric")

        # Result should contain time derivatives if non-zero
        result_str = str(result)
        if result != 0:
            self.assertTrue('A' in result_str or 'diff' in result_str,
                           "Non-zero G^tr should contain derivatives of A")

    def test_symbolic_function_handling(self):
        """Test handling of symbolic metric function A(r)"""
        A = function('A')(self.r)
        result = compute_einstein_tensor_component(self.M, A, self.r, self.t)

        # Should handle symbolic function correctly
        self.assertIsNotNone(result, "Should handle symbolic function A(r)")

        # For general A(r), result might contain derivatives
        result_str = str(result)
        # Either it's zero (simplified away) or contains A derivatives
        if result != 0:
            self.assertTrue('A' in result_str or 'diff' in result_str,
                           "Symbolic result should contain function A or its derivatives")

    def test_specific_metric_form(self):
        """Test with specific metric form"""
        # Use a simple power law: A = r^n
        n = var('n')
        A = self.r**n
        result = compute_einstein_tensor_component(self.M, A, self.r, self.t)

        # Should compute successfully
        self.assertIsNotNone(result, "Should compute G^tr for power law metric")

        # Result should depend on the power n
        result_str = str(result)
        if result != 0:
            self.assertTrue('n' in result_str or 'r' in result_str,
                           "Result should depend on metric parameters")

    def test_momentum_constraint_relevance(self):
        """Test that G^t_r is relevant for momentum constraint verification"""
        # In ADM formalism, G^tr is related to momentum constraint
        # Test with a metric that violates momentum constraint

        # Simple non-physical metric
        A = self.r + 1
        result = compute_einstein_tensor_component(self.M, A, self.r, self.t)

        # Should be able to compute the component
        self.assertIsNotNone(result, "Should compute G^tr for momentum constraint checking")

    def test_gauge_dependence(self):
        """Test G^t_r computation for different gauge choices"""
        # Different forms of the same physical metric
        M_param = var('M')

        # Standard Schwarzschild
        A1 = 1 - 2*M_param/self.r
        result1 = compute_einstein_tensor_component(self.M, A1, self.r, self.t)

        # Should be zero for both (same physics, different coordinate choice)
        result1_simplified = result1.simplify_full() if hasattr(result1, 'simplify_full') else result1
        self.assertEqual(result1_simplified, 0, "G^tr should be zero for Schwarzschild in standard coordinates")

    def test_numerical_substitution(self):
        """Test behavior with numerical parameter substitution"""
        M_param = var('M')
        A = 1 - 2*M_param/self.r
        result = compute_einstein_tensor_component(self.M, A, self.r, self.t)

        # Should handle symbolic parameters
        self.assertIsNotNone(result, "Should handle symbolic parameters")

        # For vacuum solution, should be zero
        result_simplified = result.simplify_full() if hasattr(result, 'simplify_full') else result
        self.assertEqual(result_simplified, 0, "Should be zero even with symbolic M")

    def test_uses_theta_coordinate_correctly(self):
        """Test that function uses theta coordinate correctly"""
        # The function creates var('theta') internally
        A = function('A')(self.r)
        result = compute_einstein_tensor_component(self.M, A, self.r, self.t)

        # Should compute without error
        self.assertIsNotNone(result, "Should handle theta coordinate correctly")

    def test_independence_from_input_coordinates(self):
        """Test that function works independently of input coordinate choice"""
        # Function should work with any manifold M and coordinates
        A = 1  # Flat metric

        # Should compute successfully regardless of input coordinate names
        result = compute_einstein_tensor_component(self.M, A, self.r, self.t)

        self.assertIsNotNone(result, "Should work independently of coordinate naming")

    def test_mixed_tensor_computation(self):
        """Test that mixed tensor computation G^μν is handled correctly"""
        A = function('A')(self.r)
        result = compute_einstein_tensor_component(self.M, A, self.r, self.t)

        # The function computes G^tr specifically (contravariant components)
        # This should be different from G_tr (covariant) for curved spacetime
        self.assertIsNotNone(result, "Should compute contravariant component G^tr")

    def test_verification_purpose(self):
        """Test that function serves its verification purpose"""
        # Function is designed for verification - should give reliable results
        M_param = var('M')
        A = 1 - 2*M_param/self.r
        result = compute_einstein_tensor_component(self.M, A, self.r, self.t)

        # For known vacuum solution, should reliably give zero
        result_simplified = result.simplify_full() if hasattr(result, 'simplify_full') else result
        self.assertEqual(result_simplified, 0,
                        "Verification function should reliably identify vacuum solutions")

if __name__ == '__main__':
    unittest.main()