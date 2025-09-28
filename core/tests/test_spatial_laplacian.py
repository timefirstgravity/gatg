#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Unit tests for spatial_laplacian function
Tests the spatial Laplacian operator ∇² = ∂²/∂x² + ∂²/∂y² + ∂²/∂z²
"""

import sys
import os
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sage.all import *
from core.differential_operators import spatial_laplacian

class TestSpatialLaplacian(unittest.TestCase):
    """Test spatial Laplacian operator ∇²"""

    def setUp(self):
        """Set up test fields and coordinates"""
        self.x, self.y, self.z = var('x y z')
        self.t = var('t')

        # Test fields
        self.quadratic_field = self.x**2 + self.y**2 + self.z**2
        self.polynomial_field = self.x**3 + self.y**2*self.z + self.z**3
        self.separable_field = self.x**2 * self.y**3 * self.z
        self.mixed_field = self.x*self.y + self.y*self.z + self.x*self.z

    def test_default_coordinates_quadratic_field(self):
        """Test spatial Laplacian with default coordinates on quadratic field"""
        result = spatial_laplacian(self.quadratic_field)

        # For f = x² + y² + z²:
        # ∂²f/∂x² = 2, ∂²f/∂y² = 2, ∂²f/∂z² = 2
        # ∇²f = 2 + 2 + 2 = 6
        expected = 6
        self.assertEqual(result, expected)

    def test_custom_3d_coordinates(self):
        """Test spatial Laplacian with custom 3D coordinate variables"""
        u, v, w = var('u v w')
        coords_3d = [u, v, w]

        # Field in custom coordinates
        field = u**2 + v**2 + w**2
        result = spatial_laplacian(field, coordinates=coords_3d)

        # ∇²f = ∂²f/∂u² + ∂²f/∂v² + ∂²f/∂w² = 2 + 2 + 2 = 6
        expected = 6
        self.assertEqual(result, expected)

    def test_custom_4d_coordinates_skip_time(self):
        """Test spatial Laplacian with 4D coordinates (should skip time)"""
        t, x, y, z = var('t x y z')
        coords_4d = [t, x, y, z]

        # Field that depends on all coordinates including time
        field = t**3 + x**2 + y**2 + z**2
        result = spatial_laplacian(field, coordinates=coords_4d)

        # Should only take spatial derivatives (skip time):
        # ∇²f = ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z² = 2 + 2 + 2 = 6
        # (∂²f/∂t² = 6t is ignored)
        expected = 6
        self.assertEqual(result, expected)

    def test_polynomial_field(self):
        """Test spatial Laplacian on polynomial field"""
        # f = x³ + y²z + z³
        result = spatial_laplacian(self.polynomial_field)

        # ∂²f/∂x² = ∂²(x³)/∂x² = 6x
        # ∂²f/∂y² = ∂²(y²z)/∂y² = 2z
        # ∂²f/∂z² = ∂²(z³)/∂z² = 6z
        # ∇²f = 6x + 2z + 6z = 6x + 8z
        expected = 6*self.x + 8*self.z
        self.assertEqual(result, expected)

    def test_separable_field(self):
        """Test spatial Laplacian on separable field"""
        # f = x²y³z
        result = spatial_laplacian(self.separable_field)

        # ∂²f/∂x² = ∂²(x²y³z)/∂x² = 2y³z
        # ∂²f/∂y² = ∂²(x²y³z)/∂y² = x²(6y)z = 6x²yz
        # ∂²f/∂z² = ∂²(x²y³z)/∂z² = 0 (linear in z)
        # ∇²f = 2y³z + 6x²yz
        expected = 2*self.y**3*self.z + 6*self.x**2*self.y*self.z
        self.assertEqual(result, expected)

    def test_mixed_derivative_terms(self):
        """Test spatial Laplacian on field with mixed terms"""
        # f = xy + yz + xz
        result = spatial_laplacian(self.mixed_field)

        # ∂²f/∂x² = ∂²(xy + yz + xz)/∂x² = 0
        # ∂²f/∂y² = ∂²(xy + yz + xz)/∂y² = 0
        # ∂²f/∂z² = ∂²(xy + yz + xz)/∂z² = 0
        # ∇²f = 0 (all terms are linear in each variable)
        expected = 0
        self.assertEqual(result, expected)

    def test_linearity_property(self):
        """Test linearity: ∇²(af + bg) = a∇²f + b∇²g"""
        f1 = self.x**2 + self.y**2
        f2 = self.z**2 + self.x*self.y
        a, b = var('a b')

        # Test ∇²(af + bg) = a∇²f + b∇²g
        combined_field = a*f1 + b*f2
        result_combined = spatial_laplacian(combined_field)

        result_f1 = spatial_laplacian(f1)
        result_f2 = spatial_laplacian(f2)
        expected_linear = a*result_f1 + b*result_f2

        self.assertEqual(result_combined.expand(), expected_linear.expand())

    def test_trigonometric_field(self):
        """Test spatial Laplacian on trigonometric field"""
        # f = sin(kx·r) = sin(k₁x + k₂y + k₃z)
        k1, k2, k3 = var('k1 k2 k3')
        field = sin(k1*self.x + k2*self.y + k3*self.z)

        result = spatial_laplacian(field)

        # ∂²f/∂x² = -k₁²sin(k₁x + k₂y + k₃z)
        # ∂²f/∂y² = -k₂²sin(k₁x + k₂y + k₃z)
        # ∂²f/∂z² = -k₃²sin(k₁x + k₂y + k₃z)
        # ∇²f = -(k₁² + k₂² + k₃²)sin(k₁x + k₂y + k₃z) = -|k|²f
        expected = -(k1**2 + k2**2 + k3**2) * field
        self.assertEqual(result.simplify_full(), expected.simplify_full())

    def test_exponential_field(self):
        """Test spatial Laplacian on exponential field"""
        # f = e^(ax + by + cz)
        a, b, c = var('a b c')
        field = exp(a*self.x + b*self.y + c*self.z)

        result = spatial_laplacian(field)

        # ∂²f/∂x² = a²e^(ax + by + cz) = a²f
        # ∂²f/∂y² = b²e^(ax + by + cz) = b²f
        # ∂²f/∂z² = c²e^(ax + by + cz) = c²f
        # ∇²f = (a² + b² + c²)f
        expected = (a**2 + b**2 + c**2) * field
        self.assertEqual(result.simplify_full(), expected.simplify_full())

    def test_harmonic_function(self):
        """Test spatial Laplacian on harmonic functions (∇²f = 0)"""
        # f = x² - y² is harmonic in 3D (∇²f = 0)
        harmonic_field = self.x**2 - self.y**2

        result = spatial_laplacian(harmonic_field)

        # ∂²f/∂x² = 2, ∂²f/∂y² = -2, ∂²f/∂z² = 0
        # ∇²f = 2 + (-2) + 0 = 0
        expected = 0
        self.assertEqual(result, expected)

    def test_radial_function(self):
        """Test spatial Laplacian on radial function"""
        # f = r² = x² + y² + z² (already tested, but good to be explicit)
        radial_field = self.x**2 + self.y**2 + self.z**2
        result = spatial_laplacian(radial_field)

        # For r², ∇²r² = 2 + 2 + 2 = 6 (in 3D)
        expected = 6
        self.assertEqual(result, expected)

    def test_coordinate_independence(self):
        """Test that result doesn't depend on time coordinate in 4D case"""
        t, x, y, z = var('t x y z')

        # Field with time dependence: f = t²x² + y² + z²
        field = t**2 * self.x**2 + self.y**2 + self.z**2

        # Using 4D coordinates
        coords_4d = [t, x, y, z]
        result = spatial_laplacian(field, coordinates=coords_4d)

        # Should only differentiate spatially:
        # ∂²f/∂x² = ∂²(t²x²)/∂x² = 2t²
        # ∂²f/∂y² = 2, ∂²f/∂z² = 2
        # ∇²f = 2t² + 2 + 2 = 2t² + 4
        expected = 2*t**2 + 4
        self.assertEqual(result, expected)

    def test_default_behavior_matches_explicit_coordinates(self):
        """Test that default behavior matches explicit coordinate specification"""
        field = self.x**3 + self.y**2 + self.z

        # Default coordinates
        result_default = spatial_laplacian(field)

        # Explicit 3D coordinates
        coords_explicit = [self.x, self.y, self.z]
        result_explicit = spatial_laplacian(field, coordinates=coords_explicit)

        self.assertEqual(result_default, result_explicit)

if __name__ == '__main__':
    unittest.main()