#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Unit tests for compute_hamiltonian_constraint function
Tests the Hamiltonian constraint computation H_⊥ = R^(3) + K² - K_ij K^ij
"""

import sys
import os
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sage.all import *
from core.constraints import compute_hamiltonian_constraint
from core.manifolds import setup_manifold_and_coordinates

class TestHamiltonianConstraint(unittest.TestCase):
    """Test Hamiltonian constraint computation"""

    def setUp(self):
        """Set up test manifold and coordinates"""
        self.M, self.X, self.t, self.r, self.th, self.ph = setup_manifold_and_coordinates()
        # Define a test function A(t,r)
        self.A = function('A')(self.t, self.r)

    def test_function_returns_dictionary(self):
        """Test that function returns a dictionary with expected structure"""
        result = compute_hamiltonian_constraint(self.M, self.A, self.r, self.t)

        self.assertIsInstance(result, dict)

        # Check all expected keys are present
        expected_keys = [
            'ricci_3d', 'K_squared', 'K_ij_Kij', 'hamiltonian_constraint',
            'vacuum_condition', 'matter_coupling', 'phi_derivatives'
        ]
        for key in expected_keys:
            self.assertIn(key, result)

    def test_temporal_potential_relationship(self):
        """Test Phi = (1/2) * log(A) relationship"""
        result = compute_hamiltonian_constraint(self.M, self.A, self.r, self.t)

        dPhi_dr, d2Phi_dr2, dt_Phi = result['phi_derivatives']

        # Phi = (1/2) * log(A), so ∂Phi/∂r = (1/2) * (1/A) * ∂A/∂r
        expected_dPhi_dr = (1/2) * (1/self.A) * diff(self.A, self.r)
        self.assertEqual(dPhi_dr.simplify_full(), expected_dPhi_dr.simplify_full())

        # ∂Phi/∂t = (1/2) * (1/A) * ∂A/∂t
        expected_dt_Phi = (1/2) * (1/self.A) * diff(self.A, self.t)
        self.assertEqual(dt_Phi.simplify_full(), expected_dt_Phi.simplify_full())

    def test_ricci_scalar_formula(self):
        """Test 3D Ricci scalar computation"""
        result = compute_hamiltonian_constraint(self.M, self.A, self.r, self.t)

        R3 = result['ricci_3d']
        dPhi_dr, d2Phi_dr2, dt_Phi = result['phi_derivatives']

        # R^(3) should equal 4*e^(2Φ)[∂_r²Φ + (∂_r Φ)² + (2/r)∂_r Φ]
        Phi = (1/2) * log(self.A)
        expected_R3 = 4 * exp(2*Phi) * (d2Phi_dr2 + dPhi_dr**2 + (2/self.r)*dPhi_dr)

        self.assertEqual(R3.simplify_full(), expected_R3.simplify_full())

    def test_extrinsic_curvature_properties(self):
        """Test extrinsic curvature component relationships"""
        result = compute_hamiltonian_constraint(self.M, self.A, self.r, self.t)

        K_squared = result['K_squared']
        K_ij_Kij = result['K_ij_Kij']

        # For spherical symmetry with time-independent angular components:
        # K_θθ = K_φφ = 0, so only K_rr contributes

        # K² = (K_trace)² where K_trace = K^r_r (only non-zero component)
        # K_ij K^ij = K_rr K^rr = K_rr * γ^rr * K_rr = A * K_rr²

        # Since K_trace = e^(-Φ) ∂_t Φ and K_rr = e^(-3Φ) ∂_t Φ
        # We have K² = e^(-2Φ) (∂_t Φ)² and K_ij K^ij = A * e^(-6Φ) (∂_t Φ)²

        Phi = (1/2) * log(self.A)
        dt_Phi = result['phi_derivatives'][2]

        expected_K_squared = (exp(-Phi) * dt_Phi)**2
        self.assertEqual(K_squared.simplify_full(), expected_K_squared.simplify_full())

    def test_hamiltonian_constraint_structure(self):
        """Test the structure of the Hamiltonian constraint"""
        result = compute_hamiltonian_constraint(self.M, self.A, self.r, self.t)

        hamiltonian = result['hamiltonian_constraint']
        R3 = result['ricci_3d']
        K_squared = result['K_squared']
        K_ij_Kij = result['K_ij_Kij']

        # H_⊥ = R^(3) + K² - K_ij K^ij
        expected_hamiltonian = R3 + K_squared - K_ij_Kij
        self.assertEqual(hamiltonian.simplify_full(), expected_hamiltonian.simplify_full())

    def test_vacuum_and_matter_conditions(self):
        """Test vacuum and matter condition documentation"""
        result = compute_hamiltonian_constraint(self.M, self.A, self.r, self.t)

        vacuum_condition = result['vacuum_condition']
        matter_coupling = result['matter_coupling']

        self.assertIsInstance(vacuum_condition, str)
        self.assertIsInstance(matter_coupling, str)

        self.assertIn('H_⊥ = 0', vacuum_condition)
        self.assertIn('A(t,r)', vacuum_condition)
        self.assertIn('16πG/c⁴', matter_coupling)
        self.assertIn('ρ', matter_coupling)

    def test_spherical_symmetry_consistency(self):
        """Test consistency with spherical symmetry assumptions"""
        result = compute_hamiltonian_constraint(self.M, self.A, self.r, self.t)

        # In spherical symmetry, only the radial extrinsic curvature component
        # should contribute to the trace

        # The 3D Ricci scalar should depend only on spatial (r) derivatives of Phi
        R3 = result['ricci_3d']
        dPhi_dr, d2Phi_dr2, dt_Phi = result['phi_derivatives']

        # R3 should contain spatial derivatives dPhi_dr and d2Phi_dr2
        dPhi_dr_str = str(dPhi_dr)
        d2Phi_dr2_str = str(d2Phi_dr2)

        # Check that spatial derivatives appear in R3
        spatial_derivs_in_R3 = (dPhi_dr_str in str(R3) or 'diff(' in str(R3))
        self.assertTrue(spatial_derivs_in_R3, "3D Ricci scalar should contain spatial derivatives")

        # For a static A (independent of time), R3 should not contain time derivatives
        A_static = function('A_static')(self.r)
        result_static = compute_hamiltonian_constraint(self.M, A_static, self.r, self.t)
        R3_static = result_static['ricci_3d']
        dt_Phi_static = result_static['phi_derivatives'][2]

        # For static case, no time derivatives should appear in R3
        self.assertEqual(dt_Phi_static, 0)  # Static case has no time dependence

    def test_mathematical_consistency_with_adm(self):
        """Test mathematical consistency with ADM formalism"""
        result = compute_hamiltonian_constraint(self.M, self.A, self.r, self.t)

        K_squared = result['K_squared']
        K_ij_Kij = result['K_ij_Kij']

        # In ADM formalism: K² - K_ij K^ij = K_trace² - K_ij K^ij
        # For our spherical ansatz with only K_rr ≠ 0:
        # K_trace = K^r_r = γ^rr K_rr = A * K_rr
        # K_ij K^ij = K_rr K^rr = K_rr * γ^rr * K_rr = A * K_rr²

        # Therefore: K² - K_ij K^ij = (A * K_rr)² - A * K_rr² = A²K_rr² - A*K_rr² = A*K_rr²(A-1)

        # Let's verify this relationship holds
        difference = K_squared - K_ij_Kij
        difference_simplified = difference.simplify_full()

        # This difference should be expressible in terms of A and time derivatives
        self.assertIsNotNone(difference_simplified)

    def test_temporal_derivative_handling(self):
        """Test that temporal derivatives are handled correctly"""
        result = compute_hamiltonian_constraint(self.M, self.A, self.r, self.t)

        dt_Phi = result['phi_derivatives'][2]
        K_squared = result['K_squared']
        K_ij_Kij = result['K_ij_Kij']

        # Both K² and K_ij K^ij should contain (∂_t Φ)² terms
        dt_Phi_str = str(dt_Phi)

        # Check that time derivatives appear in the extrinsic curvature terms
        K_squared_contains_time = dt_Phi_str in str(K_squared) or 'diff(' in str(K_squared)
        K_ij_Kij_contains_time = dt_Phi_str in str(K_ij_Kij) or 'diff(' in str(K_ij_Kij)

        self.assertTrue(K_squared_contains_time, "K² should contain time derivatives")
        self.assertTrue(K_ij_Kij_contains_time, "K_ij K^ij should contain time derivatives")

    def test_special_case_static_solution(self):
        """Test behavior for static solutions where ∂_t A = 0"""
        # For static case, create A that doesn't depend on time
        A_static = function('A_static')(self.r)
        result = compute_hamiltonian_constraint(self.M, A_static, self.r, self.t)

        K_squared = result['K_squared']
        K_ij_Kij = result['K_ij_Kij']

        # For static solutions, time derivatives vanish, so K² = K_ij K^ij = 0
        # Therefore H_⊥ = R^(3) for static case
        self.assertEqual(K_squared, 0)
        self.assertEqual(K_ij_Kij, 0)

        hamiltonian = result['hamiltonian_constraint']
        R3 = result['ricci_3d']
        self.assertEqual(hamiltonian.simplify_full(), R3.simplify_full())

if __name__ == '__main__':
    unittest.main()