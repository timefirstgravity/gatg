#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Unit tests for time_derivative function
Tests the time derivative operator ∂^n f/∂t^n
"""

import sys
import os
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sage.all import *
from core.differential_operators import time_derivative

class TestTimeDerivative(unittest.TestCase):
    """Test time derivative operator ∂^n f/∂t^n"""

    def setUp(self):
        """Set up test fields and coordinates"""
        self.t = var('t')
        self.x, self.y, self.z = var('x y z')

        # Test fields
        self.constant_field = 5
        self.linear_time_field = 2*self.t + 3
        self.quadratic_time_field = self.t**2 + 4*self.t + 1
        self.polynomial_field = self.t**4 + 3*self.t**3 + 2*self.t**2 + self.t + 1
        self.mixed_field = self.t**2 * self.x + self.t * self.y + self.z

    def test_first_derivative_default_order(self):
        """Test first time derivative with default order=1"""
        # f = t² + 4t + 1
        result = time_derivative(self.quadratic_time_field)

        # ∂f/∂t = 2t + 4
        expected = 2*self.t + 4
        self.assertEqual(result, expected)

    def test_first_derivative_explicit_order(self):
        """Test first time derivative with explicit order=1"""
        # f = t² + 4t + 1
        result = time_derivative(self.quadratic_time_field, order=1)

        # ∂f/∂t = 2t + 4
        expected = 2*self.t + 4
        self.assertEqual(result, expected)

    def test_zero_order_derivative(self):
        """Test zero-order derivative (should return original field)"""
        # f = t² + 4t + 1
        result = time_derivative(self.quadratic_time_field, order=0)

        # ∂⁰f/∂t⁰ = f
        expected = self.quadratic_time_field
        self.assertEqual(result, expected)

    def test_second_derivative(self):
        """Test second time derivative"""
        # f = t⁴ + 3t³ + 2t² + t + 1
        result = time_derivative(self.polynomial_field, order=2)

        # ∂²f/∂t² = 12t² + 18t + 4
        expected = 12*self.t**2 + 18*self.t + 4
        self.assertEqual(result, expected)

    def test_third_derivative(self):
        """Test third time derivative"""
        # f = t⁴ + 3t³ + 2t² + t + 1
        result = time_derivative(self.polynomial_field, order=3)

        # ∂³f/∂t³ = 24t + 18
        expected = 24*self.t + 18
        self.assertEqual(result, expected)

    def test_fourth_derivative(self):
        """Test fourth time derivative"""
        # f = t⁴ + 3t³ + 2t² + t + 1
        result = time_derivative(self.polynomial_field, order=4)

        # ∂⁴f/∂t⁴ = 24
        expected = 24
        self.assertEqual(result, expected)

    def test_higher_order_derivative_vanishes(self):
        """Test that higher order derivatives of polynomial eventually vanish"""
        # f = t⁴ + 3t³ + 2t² + t + 1 (degree 4 polynomial)
        result = time_derivative(self.polynomial_field, order=5)

        # ∂⁵f/∂t⁵ = 0 (5th derivative of degree 4 polynomial)
        expected = 0
        self.assertEqual(result, expected)

    def test_constant_field_derivative(self):
        """Test time derivative of constant field"""
        # f = 5 (constant)
        result = time_derivative(self.constant_field)

        # ∂f/∂t = 0
        expected = 0
        self.assertEqual(result, expected)

    def test_linear_time_field(self):
        """Test time derivative of linear time field"""
        # f = 2t + 3
        result = time_derivative(self.linear_time_field)

        # ∂f/∂t = 2
        expected = 2
        self.assertEqual(result, expected)

    def test_custom_time_coordinate(self):
        """Test time derivative with custom time coordinate"""
        tau = var('tau')
        field = tau**3 + 2*tau**2 + tau

        result = time_derivative(field, time_coord=tau)

        # ∂f/∂τ = 3τ² + 4τ + 1
        expected = 3*tau**2 + 4*tau + 1
        self.assertEqual(result, expected)

    def test_trigonometric_field(self):
        """Test time derivative of trigonometric field"""
        # f = sin(ωt)
        omega = var('omega')
        field = sin(omega * self.t)

        # First derivative
        result1 = time_derivative(field, order=1)
        expected1 = omega * cos(omega * self.t)
        self.assertEqual(result1, expected1)

        # Second derivative
        result2 = time_derivative(field, order=2)
        expected2 = -omega**2 * sin(omega * self.t)
        self.assertEqual(result2, expected2)

    def test_exponential_field(self):
        """Test time derivative of exponential field"""
        # f = e^(at)
        a = var('a')
        field = exp(a * self.t)

        # First derivative
        result1 = time_derivative(field, order=1)
        expected1 = a * exp(a * self.t)
        self.assertEqual(result1, expected1)

        # Second derivative
        result2 = time_derivative(field, order=2)
        expected2 = a**2 * exp(a * self.t)
        self.assertEqual(result2, expected2)

        # n-th derivative
        result_n = time_derivative(field, order=3)
        expected_n = a**3 * exp(a * self.t)
        self.assertEqual(result_n, expected_n)

    def test_mixed_spatiotemporal_field(self):
        """Test time derivative of field with both spatial and temporal dependence"""
        # f = t²x + ty + z
        result = time_derivative(self.mixed_field)

        # ∂f/∂t = 2tx + y (spatial terms unchanged)
        expected = 2*self.t*self.x + self.y
        self.assertEqual(result, expected)

    def test_logarithmic_field(self):
        """Test time derivative of logarithmic field"""
        # f = ln(t) for t > 0
        field = log(self.t)

        # First derivative
        result1 = time_derivative(field, order=1)
        expected1 = 1/self.t
        self.assertEqual(result1, expected1)

        # Second derivative
        result2 = time_derivative(field, order=2)
        expected2 = -1/self.t**2
        self.assertEqual(result2, expected2)

    def test_power_function(self):
        """Test time derivative of power function"""
        # f = t^n
        n = var('n')
        field = self.t**n

        # First derivative: ∂(t^n)/∂t = n*t^(n-1)
        result = time_derivative(field, order=1)
        expected = n * self.t**(n-1)
        self.assertEqual(result, expected)

    def test_composite_function(self):
        """Test time derivative of composite function using chain rule"""
        # f = sin(t²)
        field = sin(self.t**2)

        # First derivative: ∂sin(t²)/∂t = cos(t²) * 2t
        result = time_derivative(field, order=1)
        expected = cos(self.t**2) * 2*self.t
        self.assertEqual(result, expected)

    def test_product_rule_verification(self):
        """Test that time derivative follows product rule"""
        # f = t * sin(t), g = t, h = sin(t)
        # ∂(gh)/∂t = g'h + gh'
        g = self.t
        h = sin(self.t)
        product = g * h

        result_product = time_derivative(product)

        # Manual calculation: (t * sin(t))' = 1 * sin(t) + t * cos(t)
        expected = sin(self.t) + self.t * cos(self.t)
        self.assertEqual(result_product, expected)

    def test_iterative_application(self):
        """Test that applying time_derivative multiple times equals higher order"""
        # f = t⁴
        field = self.t**4

        # Apply first derivative twice
        first_deriv = time_derivative(field, order=1)
        second_deriv_iterative = time_derivative(first_deriv, order=1)

        # Direct second derivative
        second_deriv_direct = time_derivative(field, order=2)

        self.assertEqual(second_deriv_iterative, second_deriv_direct)

    def test_default_time_coordinate_behavior(self):
        """Test default time coordinate when none specified"""
        # Should use 't' as default time coordinate
        field = var('t')**2 + 3*var('t')

        result = time_derivative(field)  # No time_coord specified

        # ∂(t² + 3t)/∂t = 2t + 3
        expected = 2*var('t') + 3
        self.assertEqual(result, expected)

if __name__ == '__main__':
    unittest.main()