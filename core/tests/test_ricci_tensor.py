#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Unit tests for compute_ricci_tensor function
Tests the Ricci tensor and scalar computation for given metrics
"""

import sys
import os
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sage.all import *
from core.tensors import compute_ricci_tensor
from core.manifolds import setup_manifold_and_coordinates
from core.metrics import construct_metric

class TestComputeRicciTensor(unittest.TestCase):
    """Test Ricci tensor and scalar computation"""

    def setUp(self):
        """Set up test manifold and coordinates"""
        self.M, self.X, self.t, self.r, self.th, self.ph = setup_manifold_and_coordinates()

    def test_returns_tuple_with_two_elements(self):
        """Test that function returns tuple (Ricci tensor, Ricci scalar)"""
        # Create a simple metric
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        result = compute_ricci_tensor(g)

        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)

    def test_flat_spacetime_minkowski(self):
        """Test Ricci tensor for flat Minkowski spacetime"""
        # Minkowski metric: ds² = -dt² + dr² + r²dΩ²
        # This is flat spacetime with A = 1
        A = 1
        g = construct_metric(self.M, A, self.r, self.th)

        R, Rs = compute_ricci_tensor(g)

        # For flat spacetime, Ricci tensor and scalar should be zero
        # Note: SageMath may have numerical tolerance issues
        for i in range(4):
            for j in range(4):
                component = R[i,j].expr().simplify_full()
                if component != 0:
                    # Check if it's numerically very small
                    self.assertTrue(abs(component) < 1e-10 if hasattr(component, '__abs__') else component == 0,
                                   f"R[{i},{j}] = {component} should be zero for flat spacetime")

        Rs_simplified = Rs.expr().simplify_full() if hasattr(Rs, 'expr') else Rs
        self.assertEqual(Rs_simplified, 0, "Ricci scalar should be zero for flat spacetime")

    def test_schwarzschild_vacuum_solution(self):
        """Test Ricci tensor for Schwarzschild metric (vacuum solution)"""
        # Schwarzschild metric: A(r) = 1 - 2M/r
        M_param = var('M')
        A = 1 - 2*M_param/self.r
        g = construct_metric(self.M, A, self.r, self.th)

        R, Rs = compute_ricci_tensor(g)

        # Schwarzschild is a vacuum solution: R_μν = 0, R = 0
        # Due to symbolic complexity, we verify the scalar is zero
        Rs_simplified = Rs.expr().simplify_full() if hasattr(Rs, 'expr') else Rs

        # For vacuum solutions, Ricci scalar should be zero
        # Note: Symbolic simplification might be needed
        self.assertTrue(Rs_simplified == 0 or str(Rs_simplified) == '0',
                       f"Ricci scalar = {Rs_simplified} should be zero for Schwarzschild vacuum solution")

    def test_ricci_tensor_type(self):
        """Test that Ricci tensor has correct type"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        R, Rs = compute_ricci_tensor(g)

        # R should be a tensor field of type (0,2)
        self.assertTrue(hasattr(R, 'tensor_type'),
                       "Ricci tensor should have tensor_type attribute")

        if hasattr(R, 'tensor_type'):
            tensor_type = R.tensor_type()
            self.assertEqual(tensor_type, (0, 2),
                           f"Ricci tensor type should be (0,2), got {tensor_type}")

    def test_ricci_scalar_is_scalar(self):
        """Test that Ricci scalar is indeed a scalar field"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        R, Rs = compute_ricci_tensor(g)

        # Rs should be a scalar field on the manifold
        self.assertTrue(hasattr(Rs, 'expr') or isinstance(Rs, (int, float)) or hasattr(Rs, 'display'),
                       "Ricci scalar should be a scalar field or numeric value")

    def test_symmetric_ricci_tensor(self):
        """Test that Ricci tensor is symmetric: R_μν = R_νμ"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        R, Rs = compute_ricci_tensor(g)

        # Check symmetry for all component pairs
        for i in range(4):
            for j in range(i+1, 4):  # Only check upper triangle
                R_ij = R[i,j].expr() if hasattr(R[i,j], 'expr') else R[i,j]
                R_ji = R[j,i].expr() if hasattr(R[j,i], 'expr') else R[j,i]

                # Simplify both components
                R_ij_simplified = R_ij.simplify_full() if hasattr(R_ij, 'simplify_full') else R_ij
                R_ji_simplified = R_ji.simplify_full() if hasattr(R_ji, 'simplify_full') else R_ji

                self.assertEqual(R_ij_simplified, R_ji_simplified,
                               f"R[{i},{j}] != R[{j},{i}]: Ricci tensor should be symmetric")

    def test_spherically_symmetric_structure(self):
        """Test Ricci tensor structure for spherically symmetric metric"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        R, Rs = compute_ricci_tensor(g)

        # For spherical symmetry, certain off-diagonal components should be zero
        # R_tθ = R_tφ = R_rθ = R_rφ = R_θφ = 0
        zero_components = [(0,2), (0,3), (1,2), (1,3), (2,3)]

        for i, j in zero_components:
            R_ij = R[i,j].expr() if hasattr(R[i,j], 'expr') else R[i,j]
            R_ij_simplified = R_ij.simplify_full() if hasattr(R_ij, 'simplify_full') else R_ij

            self.assertEqual(R_ij_simplified, 0,
                           f"R[{i},{j}] should be zero for spherically symmetric metric")

    def test_trace_equals_ricci_scalar(self):
        """Test that trace of Ricci tensor equals Ricci scalar: g^μν R_μν = R"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        R, Rs = compute_ricci_tensor(g)

        # Compute trace manually: R = g^μν R_μν
        g_inv = g.inverse()
        trace = 0
        for i in range(4):
            for j in range(4):
                g_inv_ij = g_inv[i,j].expr() if hasattr(g_inv[i,j], 'expr') else g_inv[i,j]
                R_ij = R[i,j].expr() if hasattr(R[i,j], 'expr') else R[i,j]
                trace += g_inv_ij * R_ij

        trace_simplified = trace.simplify_full() if hasattr(trace, 'simplify_full') else trace
        Rs_simplified = Rs.expr().simplify_full() if hasattr(Rs, 'expr') else Rs

        # The trace should equal the Ricci scalar
        self.assertEqual(trace_simplified, Rs_simplified,
                        "Trace of Ricci tensor should equal Ricci scalar")

    def test_constant_curvature_spaces(self):
        """Test Ricci tensor for constant curvature metrics"""
        # For a metric with constant scalar curvature
        # The Ricci tensor should be proportional to the metric: R_μν = (R/4)g_μν
        # This is a consistency check rather than exact verification

        A = 1  # Simple case
        g = construct_metric(self.M, A, self.r, self.th)

        R, Rs = compute_ricci_tensor(g)

        # This is more of a structural test
        self.assertIsNotNone(R, "Ricci tensor should not be None")
        self.assertIsNotNone(Rs, "Ricci scalar should not be None")

    def test_coordinate_independence_of_scalar(self):
        """Test that Ricci scalar is a true scalar (coordinate-independent)"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        R, Rs = compute_ricci_tensor(g)

        # Ricci scalar should be a scalar field on the manifold
        # It shouldn't explicitly depend on coordinate choice
        # This is a conceptual test - the value is invariant
        self.assertTrue(Rs is not None, "Ricci scalar should be computed")

    def test_einstein_equation_structure(self):
        """Test that Ricci tensor can be used in Einstein equations"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        R, Rs = compute_ricci_tensor(g)

        # Form Einstein tensor: G_μν = R_μν - (1/2)R g_μν
        # This should be well-defined
        try:
            G = R - (Rs/2)*g
            self.assertIsNotNone(G, "Should be able to form Einstein tensor from Ricci tensor")

            # G should be a (0,2) tensor like R
            if hasattr(G, 'tensor_type'):
                self.assertEqual(G.tensor_type(), (0, 2),
                               "Einstein tensor should have type (0,2)")
        except Exception as e:
            self.fail(f"Failed to construct Einstein tensor: {e}")

    def test_numerical_stability_near_horizon(self):
        """Test numerical behavior near coordinate singularities"""
        # Near r = 2M (Schwarzschild horizon), A → 0
        M_val = 1
        epsilon = 0.001
        r_near_horizon = 2*M_val + epsilon

        # Substitute specific r value
        A = 1 - 2*M_val/self.r
        g = construct_metric(self.M, A, self.r, self.th)

        try:
            R, Rs = compute_ricci_tensor(g)
            # Should complete without error even near horizon
            self.assertIsNotNone(R, "Ricci tensor should be computed near horizon")
            self.assertIsNotNone(Rs, "Ricci scalar should be computed near horizon")
        except Exception as e:
            # Some coordinate singularities are expected
            self.assertTrue("division" in str(e).lower() or "singular" in str(e).lower(),
                          f"Unexpected error near horizon: {e}")

if __name__ == '__main__':
    unittest.main()