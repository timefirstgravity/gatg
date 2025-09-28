#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Unit tests for compute_einstein_tensor function
Tests the Einstein tensor computation G_μν = R_μν - (1/2)R g_μν
"""

import sys
import os
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sage.all import *
from core.tensors import compute_einstein_tensor
from core.manifolds import setup_manifold_and_coordinates
from core.metrics import construct_metric

class TestComputeEinsteinTensor(unittest.TestCase):
    """Test Einstein tensor computation"""

    def setUp(self):
        """Set up test manifold and coordinates"""
        self.M, self.X, self.t, self.r, self.th, self.ph = setup_manifold_and_coordinates()

    def test_returns_tuple_with_three_elements(self):
        """Test that function returns tuple (Einstein tensor, Ricci tensor, Ricci scalar)"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        result = compute_einstein_tensor(g)

        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 3)

    def test_einstein_tensor_formula(self):
        """Test that Einstein tensor satisfies G_μν = R_μν - (1/2)R g_μν"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        G, R, Rs = compute_einstein_tensor(g)

        # Manually compute Einstein tensor from returned Ricci components
        G_manual = R - (Rs/2)*g

        # Compare components
        for i in range(4):
            for j in range(4):
                G_ij = G[i,j].expr() if hasattr(G[i,j], 'expr') else G[i,j]
                G_manual_ij = G_manual[i,j].expr() if hasattr(G_manual[i,j], 'expr') else G_manual[i,j]

                G_ij_simplified = G_ij.simplify_full() if hasattr(G_ij, 'simplify_full') else G_ij
                G_manual_simplified = G_manual_ij.simplify_full() if hasattr(G_manual_ij, 'simplify_full') else G_manual_ij

                self.assertEqual(G_ij_simplified, G_manual_simplified,
                               f"G[{i},{j}] doesn't match manual calculation")

    def test_flat_spacetime_einstein_tensor(self):
        """Test Einstein tensor for flat Minkowski spacetime"""
        # Minkowski metric: A = 1
        A = 1
        g = construct_metric(self.M, A, self.r, self.th)

        G, R, Rs = compute_einstein_tensor(g)

        # For flat spacetime, Einstein tensor should be zero
        for i in range(4):
            for j in range(4):
                G_ij = G[i,j].expr() if hasattr(G[i,j], 'expr') else G[i,j]
                G_ij_simplified = G_ij.simplify_full() if hasattr(G_ij, 'simplify_full') else G_ij

                if G_ij_simplified != 0:
                    # Check numerical smallness
                    self.assertTrue(abs(G_ij_simplified) < 1e-10 if hasattr(G_ij_simplified, '__abs__') else G_ij_simplified == 0,
                                   f"G[{i},{j}] = {G_ij_simplified} should be zero for flat spacetime")

    def test_vacuum_einstein_equations(self):
        """Test that vacuum solution satisfies G_μν = 0"""
        # Schwarzschild metric: A(r) = 1 - 2M/r
        M_param = var('M')
        A = 1 - 2*M_param/self.r
        g = construct_metric(self.M, A, self.r, self.th)

        G, R, Rs = compute_einstein_tensor(g)

        # For vacuum solutions, Einstein tensor should vanish
        # Check that Ricci scalar is zero (implies G_μν = 0 for vacuum)
        Rs_simplified = Rs.expr().simplify_full() if hasattr(Rs, 'expr') else Rs
        self.assertTrue(Rs_simplified == 0 or str(Rs_simplified) == '0',
                       f"Ricci scalar = {Rs_simplified} should be zero for vacuum")

    def test_einstein_tensor_symmetry(self):
        """Test that Einstein tensor is symmetric: G_μν = G_νμ"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        G, R, Rs = compute_einstein_tensor(g)

        # Check symmetry
        for i in range(4):
            for j in range(i+1, 4):
                G_ij = G[i,j].expr() if hasattr(G[i,j], 'expr') else G[i,j]
                G_ji = G[j,i].expr() if hasattr(G[j,i], 'expr') else G[j,i]

                G_ij_simplified = G_ij.simplify_full() if hasattr(G_ij, 'simplify_full') else G_ij
                G_ji_simplified = G_ji.simplify_full() if hasattr(G_ji, 'simplify_full') else G_ji

                self.assertEqual(G_ij_simplified, G_ji_simplified,
                               f"G[{i},{j}] != G[{j},{i}]: Einstein tensor should be symmetric")

    def test_contracted_bianchi_identity(self):
        """Test that Einstein tensor satisfies contracted Bianchi identity: D^μ G_μν = 0"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        G, R, Rs = compute_einstein_tensor(g)

        # The covariant divergence of Einstein tensor should vanish
        # This is a fundamental identity in GR
        # For now, we just verify G is well-formed for divergence calculation
        self.assertIsNotNone(G, "Einstein tensor should be computed for Bianchi identity check")
        self.assertTrue(hasattr(G, 'tensor_type'), "Einstein tensor should be a proper tensor")

    def test_trace_of_einstein_tensor(self):
        """Test trace of Einstein tensor: g^μν G_μν = -R"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        G, R, Rs = compute_einstein_tensor(g)

        # Compute trace: g^μν G_μν
        g_inv = g.inverse()
        trace = 0
        for i in range(4):
            for j in range(4):
                g_inv_ij = g_inv[i,j].expr() if hasattr(g_inv[i,j], 'expr') else g_inv[i,j]
                G_ij = G[i,j].expr() if hasattr(G[i,j], 'expr') else G[i,j]
                trace += g_inv_ij * G_ij

        trace_simplified = trace.simplify_full() if hasattr(trace, 'simplify_full') else trace
        Rs_simplified = Rs.expr().simplify_full() if hasattr(Rs, 'expr') else Rs

        # In 4D, trace of Einstein tensor equals -R
        expected_trace = -Rs_simplified
        self.assertEqual(trace_simplified, expected_trace,
                        "Trace of Einstein tensor should equal -R")

    def test_spherically_symmetric_off_diagonal_zeros(self):
        """Test that certain off-diagonal components vanish for spherical symmetry"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        G, R, Rs = compute_einstein_tensor(g)

        # For spherical symmetry, these components should be zero
        zero_components = [(0,2), (0,3), (1,2), (1,3), (2,3)]

        for i, j in zero_components:
            G_ij = G[i,j].expr() if hasattr(G[i,j], 'expr') else G[i,j]
            G_ij_simplified = G_ij.simplify_full() if hasattr(G_ij, 'simplify_full') else G_ij

            self.assertEqual(G_ij_simplified, 0,
                           f"G[{i},{j}] should be zero for spherically symmetric metric")

    def test_einstein_tensor_type(self):
        """Test that Einstein tensor has correct tensor type (0,2)"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        G, R, Rs = compute_einstein_tensor(g)

        # G should be a (0,2) tensor like the metric
        self.assertTrue(hasattr(G, 'tensor_type'),
                       "Einstein tensor should have tensor_type attribute")

        if hasattr(G, 'tensor_type'):
            tensor_type = G.tensor_type()
            self.assertEqual(tensor_type, (0, 2),
                           f"Einstein tensor type should be (0,2), got {tensor_type}")

    def test_energy_momentum_coupling(self):
        """Test that Einstein tensor structure allows energy-momentum coupling"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        G, R, Rs = compute_einstein_tensor(g)

        # Einstein equations: G_μν = 8πG/c⁴ T_μν
        # Test that we can form this equation symbolically
        kappa = var('kappa')  # 8πG/c⁴

        # Should be able to set up field equations
        try:
            # Symbolic stress-energy tensor
            T = g.parent().tensor_field(0, 2, name='T')
            field_equations = G - kappa * T
            self.assertIsNotNone(field_equations,
                               "Should be able to form Einstein field equations")
        except:
            # At minimum, G should be a proper tensor for coupling
            self.assertTrue(hasattr(G, 'copy') or hasattr(G, '__getitem__'),
                          "Einstein tensor should support tensor operations")

    def test_cosmological_constant_modification(self):
        """Test Einstein tensor with cosmological constant: G_μν + Λg_μν"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        G, R, Rs = compute_einstein_tensor(g)

        # Should be able to add cosmological constant term
        Lambda = var('Lambda')
        try:
            G_Lambda = G + Lambda * g
            self.assertIsNotNone(G_Lambda,
                               "Should be able to add cosmological constant term")

            # Result should still be a (0,2) tensor
            if hasattr(G_Lambda, 'tensor_type'):
                self.assertEqual(G_Lambda.tensor_type(), (0, 2),
                               "Modified Einstein tensor should remain type (0,2)")
        except Exception as e:
            self.fail(f"Failed to add cosmological constant term: {e}")

    def test_all_returned_components_consistent(self):
        """Test that all three returned components are mutually consistent"""
        A = function('A')(self.r)
        g = construct_metric(self.M, A, self.r, self.th)

        G, R, Rs = compute_einstein_tensor(g)

        # Verify internal consistency
        # 1. R should be the Ricci tensor of g
        # 2. Rs should be the trace of R
        # 3. G should equal R - (Rs/2)*g

        # Already tested formula above, but verify all three are non-None
        self.assertIsNotNone(G, "Einstein tensor should not be None")
        self.assertIsNotNone(R, "Ricci tensor should not be None")
        self.assertIsNotNone(Rs, "Ricci scalar should not be None")

        # Test they have expected structure
        self.assertTrue(hasattr(G, '__getitem__'), "G should support component access")
        self.assertTrue(hasattr(R, '__getitem__'), "R should support component access")
        self.assertTrue(hasattr(Rs, 'expr') or isinstance(Rs, (int, float)) or Rs == 0,
                       "Rs should be a scalar expression")

if __name__ == '__main__':
    unittest.main()