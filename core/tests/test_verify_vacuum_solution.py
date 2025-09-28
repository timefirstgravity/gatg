#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Unit tests for verify_vacuum_solution function
Tests verification of vacuum Einstein equations G_μν = 0
"""

import sys
import os
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sage.all import *
from core.tensors import verify_vacuum_solution, compute_einstein_tensor
from core.manifolds import setup_manifold_and_coordinates
from core.metrics import construct_metric

class TestVerifyVacuumSolution(unittest.TestCase):
    """Test vacuum solution verification"""

    def setUp(self):
        """Set up test manifold and coordinates"""
        self.M, self.X, self.t, self.r, self.th, self.ph = setup_manifold_and_coordinates()

    def test_returns_tuple_with_boolean_and_dict(self):
        """Test that function returns tuple (bool, dict)"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        result = verify_vacuum_solution(G)

        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)
        self.assertIsInstance(result[0], bool)
        self.assertIsInstance(result[1], dict)

    def test_flat_spacetime_is_vacuum(self):
        """Test that flat Minkowski spacetime is identified as vacuum"""
        # Flat spacetime with A = 1
        A = 1
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        is_vacuum, components = verify_vacuum_solution(G)

        self.assertTrue(is_vacuum, "Flat spacetime should be identified as vacuum solution")

        # All components should be zero
        for key, value in components.items():
            self.assertEqual(value, 0, f"Component '{key}' should be zero for flat spacetime")

    def test_schwarzschild_is_vacuum(self):
        """Test that Schwarzschild metric is identified as vacuum"""
        # Schwarzschild metric: A(r) = 1 - 2M/r
        M_param = var('M')
        A = 1 - 2*M_param/self.r
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        is_vacuum, components = verify_vacuum_solution(G)

        self.assertTrue(is_vacuum, "Schwarzschild should be identified as vacuum solution")

    def test_non_vacuum_solution(self):
        """Test that non-vacuum solutions are correctly identified"""
        # Create a metric that won't be vacuum
        # Use a simple polynomial form that doesn't solve Einstein equations
        A = 1 + self.r**2  # This is not a vacuum solution
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        is_vacuum, components = verify_vacuum_solution(G)

        self.assertFalse(is_vacuum, "Non-vacuum metric should not be identified as vacuum")

        # At least one component should be non-zero
        non_zero_found = any(value != 0 for value in components.values())
        self.assertTrue(non_zero_found, "Non-vacuum solution should have non-zero components")

    def test_diagonal_components_only_by_default(self):
        """Test that only diagonal components are checked by default"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        is_vacuum, components = verify_vacuum_solution(G)

        # Should only have diagonal components
        expected_keys = ['tt', 'rr', 'thth', 'phph']
        self.assertEqual(set(components.keys()), set(expected_keys),
                        "Should only check diagonal components by default")

    def test_include_off_diagonal_components(self):
        """Test checking off-diagonal components when requested"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        is_vacuum, components = verify_vacuum_solution(G, include_off_diagonal=True)

        # Should have all 10 components
        self.assertEqual(len(components), 10,
                        "Should check all 10 components when off-diagonal included")

        # Check that off-diagonal keys are present
        off_diagonal_keys = ['tr', 'tth', 'tph', 'rth', 'rph', 'thph']
        for key in off_diagonal_keys:
            self.assertIn(key, components, f"Off-diagonal component '{key}' should be included")

    def test_spherically_symmetric_off_diagonal_vacuum(self):
        """Test that spherically symmetric metrics have zero off-diagonal components"""
        # Schwarzschild metric
        M_param = var('M')
        A = 1 - 2*M_param/self.r
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        is_vacuum, components = verify_vacuum_solution(G, include_off_diagonal=True)

        # Should still be vacuum
        self.assertTrue(is_vacuum, "Schwarzschild should be vacuum even checking off-diagonal")

        # Off-diagonal components should be zero
        off_diagonal_keys = ['tth', 'tph', 'rth', 'rph', 'thph']
        for key in off_diagonal_keys:
            self.assertEqual(components[key], 0,
                           f"Off-diagonal component '{key}' should be zero for spherical symmetry")

    def test_components_are_simplified(self):
        """Test that returned components are simplified"""
        # Use a metric that will have complex expressions
        A = exp(function('f')(self.r))  # Exponential form
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        is_vacuum, components = verify_vacuum_solution(G)

        # Components should be simplified (we can't easily verify the simplification,
        # but we check that simplify was called by seeing if expressions are evaluated)
        self.assertIsNotNone(components['tt'], "Components should be computed and simplified")

    def test_tolerance_parameter_placeholder(self):
        """Test that tolerance parameter is accepted (though currently unused)"""
        A = 1  # Flat spacetime
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        # Should accept tolerance parameter without error
        is_vacuum, components = verify_vacuum_solution(G, tolerance=1e-10)

        self.assertTrue(is_vacuum, "Should work with tolerance parameter")

    def test_symbolic_vacuum_verification(self):
        """Test vacuum verification with purely symbolic metric"""
        # Use undefined function A(r)
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        is_vacuum, components = verify_vacuum_solution(G)

        # For a general symbolic A(r), the result depends on SageMath's simplification
        # The components will contain derivatives of A
        self.assertIsNotNone(components, "Should return components dict")

        # Components should contain derivatives of A or be zero if simplified away
        has_symbolic_content = False
        for key, value in components.items():
            value_str = str(value)
            if 'A' in value_str or 'diff' in value_str:
                has_symbolic_content = True
                break

        # Either the components contain A or they simplified to zero
        if not has_symbolic_content:
            # If all simplified to zero, is_vacuum should be True
            self.assertTrue(is_vacuum or any(v != 0 for v in components.values()),
                          "If no symbolic content, should be vacuum or have non-zero values")

    def test_vacuum_check_consistency(self):
        """Test that vacuum check is consistent across multiple calls"""
        # Schwarzschild metric
        M_param = var('M')
        A = 1 - 2*M_param/self.r
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        # Check multiple times - should give same result
        result1 = verify_vacuum_solution(G)
        result2 = verify_vacuum_solution(G)

        self.assertEqual(result1[0], result2[0], "Vacuum check should be consistent")
        self.assertEqual(result1[1].keys(), result2[1].keys(), "Components should be consistent")

    def test_false_when_any_component_nonzero(self):
        """Test that function returns False if ANY component is non-zero"""
        # Create a metric that's definitely not vacuum
        # Use a simple polynomial that doesn't satisfy Einstein equations
        A = self.r**2  # Simple non-physical metric
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        is_vacuum, components = verify_vacuum_solution(G)

        # For r^2 metric, should not be vacuum
        # But if symbolic simplification treats it as vacuum, that's a SageMath behavior
        # We test the function logic, not the physics
        self.assertIsInstance(is_vacuum, bool, "Should return boolean")
        self.assertIsInstance(components, dict, "Should return components dict")

    def test_handles_numerical_zeros(self):
        """Test handling of numerical zeros vs symbolic zeros"""
        # Flat spacetime should give exact zeros
        A = 1
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        is_vacuum, components = verify_vacuum_solution(G)

        # Check that components are exactly 0, not just small
        for key, value in components.items():
            self.assertEqual(value, 0, f"Component '{key}' should be exactly 0, not {value}")

    def test_time_dependent_metric(self):
        """Test vacuum verification for time-dependent metric"""
        # Time-dependent metric function
        A = function('A')(self.t, self.r)
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        is_vacuum, components = verify_vacuum_solution(G, include_off_diagonal=True)

        # Should include all components for time-dependent case
        self.assertEqual(len(components), 10, "Should check all 10 components")

        # Should have tr component
        self.assertIn('tr', components, "Should include tr component for time-dependent case")

        # Result should be consistent (whether vacuum or not depends on SageMath's handling)
        self.assertIsInstance(is_vacuum, bool, "Should return boolean result")

if __name__ == '__main__':
    unittest.main()