#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Unit tests for dalembertian function
Tests the D'Alembertian wave operator □ = η^μν ∂_μ ∂_ν
"""

import sys
import os
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sage.all import *
from core.differential_operators import dalembertian

class TestDalembertian(unittest.TestCase):
    """Test D'Alembertian wave operator"""

    def setUp(self):
        """Set up test fields and coordinates"""
        self.t, self.x, self.y, self.z = var('t x y z')
        self.coords = [self.t, self.x, self.y, self.z]

        # Test fields
        self.scalar_field = self.t**2 + self.x**2 + self.y**2 + self.z**2
        self.wave_field = sin(self.t - self.x)  # Simple wave solution
        self.polynomial_field = self.t**3 + self.x*self.y + self.z**2

    def test_lorentzian_signature_default_coordinates(self):
        """Test D'Alembertian with Lorentzian signature and default coordinates"""
        result = dalembertian(self.scalar_field, signature='lorentzian')

        # For field f = t² + x² + y² + z², second derivatives are all 2
        # □f = -∂²f/∂t² + ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z² = -2 + 2 + 2 + 2 = 4
        expected = 4
        self.assertEqual(result, expected)

    def test_euclidean_signature_default_coordinates(self):
        """Test D'Alembertian with Euclidean signature and default coordinates"""
        result = dalembertian(self.scalar_field, signature='euclidean')

        # For field f = t² + x² + y² + z², all second derivatives are 2
        # □f = ∂²f/∂t² + ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z² = 2 + 2 + 2 + 2 = 8
        expected = 8
        self.assertEqual(result, expected)

    def test_custom_coordinates(self):
        """Test D'Alembertian with custom coordinate variables"""
        u, v, w, s = var('u v w s')
        custom_coords = [u, v, w, s]

        # Field in custom coordinates
        field = u**2 + v**2 + w**2 + s**2

        result = dalembertian(field, coordinates=custom_coords, signature='lorentzian')

        # □f = -∂²f/∂u² + ∂²f/∂v² + ∂²f/∂w² + ∂²f/∂s² = -2 + 2 + 2 + 2 = 4
        expected = 4
        self.assertEqual(result, expected)

    def test_wave_equation_solution(self):
        """Test D'Alembertian on a known wave equation solution"""
        # f(t,x) = sin(t - x) is a solution to the 1D wave equation □f = 0
        # In 4D: f(t,x,y,z) = sin(t - x) (independent of y,z)
        wave_field = sin(self.t - self.x)

        result = dalembertian(wave_field, signature='lorentzian')

        # ∂²f/∂t² = -sin(t-x), ∂²f/∂x² = -sin(t-x), ∂²f/∂y² = 0, ∂²f/∂z² = 0
        # □f = -(-sin(t-x)) + (-sin(t-x)) + 0 + 0 = sin(t-x) - sin(t-x) = 0
        expected = 0
        self.assertEqual(result.simplify_full(), expected)

    def test_polynomial_field(self):
        """Test D'Alembertian on polynomial field"""
        # f = t³ + xy + z²
        field = self.t**3 + self.x*self.y + self.z**2

        result = dalembertian(field, signature='lorentzian')

        # ∂²f/∂t² = 6t, ∂²f/∂x² = 0, ∂²f/∂y² = 0, ∂²f/∂z² = 2
        # □f = -6t + 0 + 0 + 2 = 2 - 6t
        expected = 2 - 6*self.t
        self.assertEqual(result, expected)

    def test_linearity_property(self):
        """Test linearity: □(af + bg) = a□f + b□g"""
        f1 = self.t**2 + self.x**2
        f2 = self.y**2 + self.z**2
        a, b = var('a b')

        # Test □(af + bg) = a□f + b□g
        combined_field = a*f1 + b*f2
        result_combined = dalembertian(combined_field, signature='lorentzian')

        result_f1 = dalembertian(f1, signature='lorentzian')
        result_f2 = dalembertian(f2, signature='lorentzian')
        expected_linear = a*result_f1 + b*result_f2

        self.assertEqual(result_combined.expand(), expected_linear.expand())

    def test_signature_difference(self):
        """Test difference between Lorentzian and Euclidean signatures"""
        field = self.t**2

        lorentzian_result = dalembertian(field, signature='lorentzian')
        euclidean_result = dalembertian(field, signature='euclidean')

        # For f = t², ∂²f/∂t² = 2
        # Lorentzian: □f = -2, Euclidean: □f = +2
        self.assertEqual(lorentzian_result, -2)
        self.assertEqual(euclidean_result, 2)

    def test_invalid_signature_raises_error(self):
        """Test that invalid signature raises ValueError"""
        with self.assertRaises(ValueError) as context:
            dalembertian(self.scalar_field, signature='invalid')

        self.assertIn("Unknown signature", str(context.exception))
        self.assertIn("invalid", str(context.exception))

    def test_coordinate_truncation(self):
        """Test that only first 4 coordinates are used"""
        # Provide more than 4 coordinates
        extra_coords = [self.t, self.x, self.y, self.z, var('w'), var('v')]

        field = self.t**2 + self.x**2 + self.y**2 + self.z**2
        result = dalembertian(field, coordinates=extra_coords, signature='lorentzian')

        # Should still work with first 4 coordinates
        expected = 4  # Same as default case
        self.assertEqual(result, expected)

    def test_mixed_derivative_terms(self):
        """Test field with mixed derivative terms"""
        # Field with cross terms: f = txy
        field = self.t * self.x * self.y

        result = dalembertian(field, signature='lorentzian')

        # ∂²f/∂t² = 0, ∂²f/∂x² = 0, ∂²f/∂y² = 0, ∂²f/∂z² = 0
        # □f = 0
        expected = 0
        self.assertEqual(result, expected)

    def test_exponential_field(self):
        """Test D'Alembertian on exponential field"""
        # f = e^(at + bx + cy + dz)
        a, b, c, d = var('a b c d')
        field = exp(a*self.t + b*self.x + c*self.y + d*self.z)

        result = dalembertian(field, signature='lorentzian')

        # ∂²f/∂t² = a²f, ∂²f/∂x² = b²f, ∂²f/∂y² = c²f, ∂²f/∂z² = d²f
        # □f = (-a² + b² + c² + d²)f
        expected = (-a**2 + b**2 + c**2 + d**2) * field
        self.assertEqual(result.simplify_full(), expected.simplify_full())

    def test_trigonometric_field(self):
        """Test D'Alembertian on trigonometric field"""
        # f = cos(ωt + kx)
        omega, k = var('omega k')
        field = cos(omega*self.t + k*self.x)

        result = dalembertian(field, signature='lorentzian')

        # ∂²f/∂t² = -ω²cos(ωt+kx), ∂²f/∂x² = -k²cos(ωt+kx)
        # □f = -(-ω²) + (-k²) + 0 + 0 = (ω² - k²)cos(ωt+kx)
        expected = (omega**2 - k**2) * field
        self.assertEqual(result.simplify_full(), expected.simplify_full())

    def test_default_behavior_no_coordinates(self):
        """Test default behavior when no coordinates provided"""
        # Should use default t, x, y, z variables
        field = var('t')**2 + var('x')**2 + var('y')**2 + var('z')**2

        result = dalembertian(field)  # Default signature is lorentzian

        # Same as lorentzian test
        expected = 4
        self.assertEqual(result, expected)

if __name__ == '__main__':
    unittest.main()