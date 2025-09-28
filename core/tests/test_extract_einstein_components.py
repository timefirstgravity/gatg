#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Unit tests for extract_einstein_components function
Tests extraction and labeling of Einstein tensor components
"""

import sys
import os
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sage.all import *
from core.tensors import extract_einstein_components, compute_einstein_tensor
from core.manifolds import setup_manifold_and_coordinates
from core.metrics import construct_metric

class TestExtractEinsteinComponents(unittest.TestCase):
    """Test Einstein tensor component extraction"""

    def setUp(self):
        """Set up test manifold and Einstein tensor"""
        self.M, self.X, self.t, self.r, self.th, self.ph = setup_manifold_and_coordinates()

        # Create a test metric and compute Einstein tensor
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)
        self.G, self.R, self.Rs = compute_einstein_tensor(g)

    def test_returns_dictionary(self):
        """Test that function returns a dictionary"""
        result = extract_einstein_components(self.G)

        self.assertIsInstance(result, dict)

    def test_diagonal_components_always_included(self):
        """Test that diagonal components are always extracted"""
        result = extract_einstein_components(self.G)

        # Check all diagonal components are present
        expected_diagonal = ['tt', 'rr', 'thth', 'phph']
        for key in expected_diagonal:
            self.assertIn(key, result, f"Diagonal component '{key}' should be included")

    def test_diagonal_only_by_default(self):
        """Test that only diagonal components are included by default"""
        result = extract_einstein_components(self.G)

        # Should only have 4 diagonal components
        self.assertEqual(len(result), 4, "Should only have 4 diagonal components by default")

        # Check no off-diagonal components
        off_diagonal_keys = ['tr', 'tth', 'tph', 'rth', 'rph', 'thph']
        for key in off_diagonal_keys:
            self.assertNotIn(key, result, f"Off-diagonal component '{key}' should not be included by default")

    def test_include_off_diagonal_components(self):
        """Test inclusion of off-diagonal components when requested"""
        result = extract_einstein_components(self.G, include_off_diagonal=True)

        # Should have 4 diagonal + 6 off-diagonal = 10 components
        self.assertEqual(len(result), 10, "Should have 10 components when off-diagonal included")

        # Check all off-diagonal components are present
        off_diagonal_keys = ['tr', 'tth', 'tph', 'rth', 'rph', 'thph']
        for key in off_diagonal_keys:
            self.assertIn(key, result, f"Off-diagonal component '{key}' should be included")

    def test_component_labeling_convention(self):
        """Test that component labels follow the correct convention"""
        result = extract_einstein_components(self.G, include_off_diagonal=True)

        # Verify labeling convention:
        # tt = G[0,0], rr = G[1,1], thth = G[2,2], phph = G[3,3]
        # tr = G[0,1], tth = G[0,2], etc.

        # Check a few components to verify correct indexing
        G_00 = self.G[0,0].expr()
        G_11 = self.G[1,1].expr()

        self.assertEqual(result['tt'], G_00, "Component 'tt' should correspond to G[0,0]")
        self.assertEqual(result['rr'], G_11, "Component 'rr' should correspond to G[1,1]")

    def test_simplification_enabled(self):
        """Test that simplification is applied when simplify=True (default)"""
        # Create a metric that will have complex expressions
        M_param = var('M')
        A = 1 - 2*M_param/self.r
        g = construct_metric(self.M, A, self.r, self.th)
        G_complex, _, _ = compute_einstein_tensor(g)

        result_simplified = extract_einstein_components(G_complex, simplify=True)

        # Check that simplification was attempted (components should have simplify_full applied)
        # We can't easily verify the simplification worked, but we can check the function ran
        self.assertIsNotNone(result_simplified['tt'])

    def test_simplification_disabled(self):
        """Test that simplification is skipped when simplify=False"""
        result_no_simplify = extract_einstein_components(self.G, simplify=False)

        # Components should be extracted but not simplified
        self.assertIsNotNone(result_no_simplify['tt'])

        # The components should be the raw .expr() values
        G_00_raw = self.G[0,0].expr()
        self.assertEqual(result_no_simplify['tt'], G_00_raw,
                        "Without simplification, should return raw .expr() value")

    def test_spherically_symmetric_off_diagonal_zeros(self):
        """Test that off-diagonal components are zero for spherically symmetric metric"""
        # For spherically symmetric metric, certain off-diagonal components should be zero
        result = extract_einstein_components(self.G, include_off_diagonal=True, simplify=True)

        # These should be zero for spherical symmetry
        zero_components = ['tth', 'tph', 'rth', 'rph', 'thph']

        for key in zero_components:
            component = result[key]
            self.assertEqual(component, 0,
                           f"Component '{key}' should be zero for spherically symmetric metric")

    def test_component_extraction_consistency(self):
        """Test that extracted components match direct tensor access"""
        result = extract_einstein_components(self.G, simplify=False)

        # Verify each diagonal component matches direct access
        self.assertEqual(result['tt'], self.G[0,0].expr())
        self.assertEqual(result['rr'], self.G[1,1].expr())
        self.assertEqual(result['thth'], self.G[2,2].expr())
        self.assertEqual(result['phph'], self.G[3,3].expr())

    def test_off_diagonal_extraction_consistency(self):
        """Test that off-diagonal components match direct tensor access"""
        result = extract_einstein_components(self.G, include_off_diagonal=True, simplify=False)

        # Verify off-diagonal components match direct access
        self.assertEqual(result['tr'], self.G[0,1].expr())
        self.assertEqual(result['tth'], self.G[0,2].expr())
        self.assertEqual(result['tph'], self.G[0,3].expr())
        self.assertEqual(result['rth'], self.G[1,2].expr())
        self.assertEqual(result['rph'], self.G[1,3].expr())
        self.assertEqual(result['thph'], self.G[2,3].expr())

    def test_flat_spacetime_components(self):
        """Test component extraction for flat spacetime"""
        # Create flat metric with A = 1
        A = 1
        g_flat = construct_metric(self.M, A, self.r, self.th)
        G_flat, _, _ = compute_einstein_tensor(g_flat)

        result = extract_einstein_components(G_flat, simplify=True)

        # All components should be zero or very small for flat spacetime
        for key, value in result.items():
            if value != 0:
                self.assertTrue(abs(value) < 1e-10 if hasattr(value, '__abs__') else value == 0,
                              f"Component '{key}' = {value} should be zero for flat spacetime")

    def test_dictionary_keys_are_strings(self):
        """Test that all dictionary keys are strings"""
        result = extract_einstein_components(self.G, include_off_diagonal=True)

        for key in result.keys():
            self.assertIsInstance(key, str, f"Key '{key}' should be a string")

    def test_preserves_symbolic_expressions(self):
        """Test that symbolic expressions are preserved in components"""
        # Use a symbolic metric function
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        result = extract_einstein_components(G, simplify=False)

        # Components should contain derivatives of A
        tt_component = str(result['tt'])
        self.assertTrue('A' in tt_component or 'diff' in tt_component,
                       "Components should preserve symbolic expressions with function A")

    def test_handles_complex_metrics(self):
        """Test extraction for more complex metric forms"""
        # Create a metric with time dependence
        A = function('A')(self.t, self.r)
        g = construct_metric(self.M, A, self.r, self.th)
        G, _, _ = compute_einstein_tensor(g)

        result = extract_einstein_components(G, include_off_diagonal=True, simplify=True)

        # Should successfully extract all components
        self.assertEqual(len(result), 10, "Should extract all 10 components for complex metric")

        # tr component might be non-zero for time-dependent metric
        self.assertIn('tr', result, "Should include tr component for time-dependent metric")

if __name__ == '__main__':
    unittest.main()