#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Unit tests for gradient_3d function
Tests the 3D gradient operator ∇f = (∂f/∂x, ∂f/∂y, ∂f/∂z)
"""

import sys
import os
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sage.all import *
from core.differential_operators import gradient_3d

class TestGradient3D(unittest.TestCase):
    """Test 3D gradient operator ∇f"""

    def setUp(self):
        """Set up test fields and coordinates"""
        self.x, self.y, self.z = var('x y z')
        self.t = var('t')

        # Test fields
        self.linear_field = 2*self.x + 3*self.y + 4*self.z
        self.quadratic_field = self.x**2 + self.y**2 + self.z**2
        self.polynomial_field = self.x**3 + self.y**2*self.z + self.z**3
        self.separable_field = self.x**2 * self.y * self.z**3

    def test_returns_list_of_three_elements(self):
        """Test that gradient returns a list with exactly 3 elements"""
        result = gradient_3d(self.quadratic_field)

        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 3)

    def test_linear_field_default_coordinates(self):
        """Test gradient of linear field with default coordinates"""
        # f = 2x + 3y + 4z
        result = gradient_3d(self.linear_field)

        # ∇f = (∂f/∂x, ∂f/∂y, ∂f/∂z) = (2, 3, 4)
        expected = [2, 3, 4]
        self.assertEqual(result, expected)

    def test_quadratic_field_default_coordinates(self):
        """Test gradient of quadratic field with default coordinates"""
        # f = x² + y² + z²
        result = gradient_3d(self.quadratic_field)

        # ∇f = (∂f/∂x, ∂f/∂y, ∂f/∂z) = (2x, 2y, 2z)
        expected = [2*self.x, 2*self.y, 2*self.z]
        self.assertEqual(result, expected)

    def test_custom_3d_coordinates(self):
        """Test gradient with custom 3D coordinate variables"""
        u, v, w = var('u v w')
        coords_3d = [u, v, w]

        # Field in custom coordinates: f = u² + v² + w²
        field = u**2 + v**2 + w**2
        result = gradient_3d(field, coordinates=coords_3d)

        # ∇f = (∂f/∂u, ∂f/∂v, ∂f/∂w) = (2u, 2v, 2w)
        expected = [2*u, 2*v, 2*w]
        self.assertEqual(result, expected)

    def test_custom_4d_coordinates_skip_time(self):
        """Test gradient with 4D coordinates (should skip time)"""
        t, x, y, z = var('t x y z')
        coords_4d = [t, x, y, z]

        # Field depending on all coordinates including time
        field = t**2 + x**2 + y**2 + z**2
        result = gradient_3d(field, coordinates=coords_4d)

        # Should only compute spatial derivatives:
        # ∇f = (∂f/∂x, ∂f/∂y, ∂f/∂z) = (2x, 2y, 2z)
        # (∂f/∂t = 2t is not included)
        expected = [2*x, 2*y, 2*z]
        self.assertEqual(result, expected)

    def test_polynomial_field(self):
        """Test gradient on polynomial field"""
        # f = x³ + y²z + z³
        result = gradient_3d(self.polynomial_field)

        # ∂f/∂x = 3x²
        # ∂f/∂y = 2yz
        # ∂f/∂z = y² + 3z²
        expected = [3*self.x**2, 2*self.y*self.z, self.y**2 + 3*self.z**2]
        self.assertEqual(result, expected)

    def test_separable_field(self):
        """Test gradient on separable field"""
        # f = x²yz³
        result = gradient_3d(self.separable_field)

        # ∂f/∂x = 2xyz³
        # ∂f/∂y = x²z³
        # ∂f/∂z = 3x²yz²
        expected = [
            2*self.x*self.y*self.z**3,
            self.x**2*self.z**3,
            3*self.x**2*self.y*self.z**2
        ]
        self.assertEqual(result, expected)

    def test_constant_field(self):
        """Test gradient of constant field"""
        constant = var('c')
        result = gradient_3d(constant)

        # ∇c = (0, 0, 0) for any constant
        expected = [0, 0, 0]
        self.assertEqual(result, expected)

    def test_trigonometric_field(self):
        """Test gradient on trigonometric field"""
        # f = sin(x) + cos(y) + sin(z)
        field = sin(self.x) + cos(self.y) + sin(self.z)
        result = gradient_3d(field)

        # ∂f/∂x = cos(x)
        # ∂f/∂y = -sin(y)
        # ∂f/∂z = cos(z)
        expected = [cos(self.x), -sin(self.y), cos(self.z)]
        self.assertEqual(result, expected)

    def test_exponential_field(self):
        """Test gradient on exponential field"""
        # f = e^(x + y + z)
        field = exp(self.x + self.y + self.z)
        result = gradient_3d(field)

        # ∂f/∂x = e^(x+y+z) = f
        # ∂f/∂y = e^(x+y+z) = f
        # ∂f/∂z = e^(x+y+z) = f
        expected = [field, field, field]
        self.assertEqual(result, expected)

    def test_logarithmic_field(self):
        """Test gradient on logarithmic field"""
        # f = ln(xyz) = ln(x) + ln(y) + ln(z)
        field = log(self.x) + log(self.y) + log(self.z)
        result = gradient_3d(field)

        # ∂f/∂x = 1/x
        # ∂f/∂y = 1/y
        # ∂f/∂z = 1/z
        expected = [1/self.x, 1/self.y, 1/self.z]
        self.assertEqual(result, expected)

    def test_mixed_product_field(self):
        """Test gradient on field with mixed products"""
        # f = xy + yz + xz
        field = self.x*self.y + self.y*self.z + self.x*self.z
        result = gradient_3d(field)

        # ∂f/∂x = y + z
        # ∂f/∂y = x + z
        # ∂f/∂z = y + x
        expected = [self.y + self.z, self.x + self.z, self.y + self.x]
        self.assertEqual(result, expected)

    def test_radial_function(self):
        """Test gradient of radial function r = √(x² + y² + z²)"""
        # f = √(x² + y² + z²)
        field = sqrt(self.x**2 + self.y**2 + self.z**2)
        result = gradient_3d(field)

        # ∂f/∂x = x/√(x² + y² + z²) = x/r
        # ∂f/∂y = y/√(x² + y² + z²) = y/r
        # ∂f/∂z = z/√(x² + y² + z²) = z/r
        r = field
        expected = [self.x/r, self.y/r, self.z/r]

        # Simplify both sides for comparison
        result_simplified = [expr.simplify_full() for expr in result]
        expected_simplified = [expr.simplify_full() for expr in expected]
        self.assertEqual(result_simplified, expected_simplified)

    def test_linearity_property(self):
        """Test linearity: ∇(af + bg) = a∇f + b∇g"""
        f1 = self.x**2 + self.y
        f2 = self.z**2 + self.x
        a, b = var('a b')

        # Test ∇(af + bg) = a∇f + b∇g
        combined_field = a*f1 + b*f2
        result_combined = gradient_3d(combined_field)

        grad_f1 = gradient_3d(f1)
        grad_f2 = gradient_3d(f2)
        expected_linear = [a*grad_f1[i] + b*grad_f2[i] for i in range(3)]

        # Expand both sides for comparison
        result_expanded = [expr.expand() for expr in result_combined]
        expected_expanded = [expr.expand() for expr in expected_linear]
        self.assertEqual(result_expanded, expected_expanded)

    def test_time_independence_in_4d_case(self):
        """Test that time coordinate doesn't affect spatial gradient"""
        t, x, y, z = var('t x y z')

        # Field with time dependence: f = t³ + x² + y² + z²
        field = t**3 + self.x**2 + self.y**2 + self.z**2

        # Using 4D coordinates should ignore time derivatives
        coords_4d = [t, x, y, z]
        result = gradient_3d(field, coordinates=coords_4d)

        # Should only compute spatial derivatives: (2x, 2y, 2z)
        expected = [2*x, 2*y, 2*z]
        self.assertEqual(result, expected)

    def test_default_behavior_matches_explicit_coordinates(self):
        """Test that default behavior matches explicit coordinate specification"""
        field = self.x**3 + self.y**2 + self.z

        # Default coordinates
        result_default = gradient_3d(field)

        # Explicit 3D coordinates
        coords_explicit = [self.x, self.y, self.z]
        result_explicit = gradient_3d(field, coordinates=coords_explicit)

        self.assertEqual(result_default, result_explicit)

if __name__ == '__main__':
    unittest.main()