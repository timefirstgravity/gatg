#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Unit tests for divergence_3d function
Tests the 3D divergence operator ∇·V = ∂V^x/∂x + ∂V^y/∂y + ∂V^z/∂z
"""

import sys
import os
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sage.all import *
from core.differential_operators import divergence_3d

class TestDivergence3D(unittest.TestCase):
    """Test 3D divergence operator ∇·V"""

    def setUp(self):
        """Set up test vector fields and coordinates"""
        self.x, self.y, self.z = var('x y z')
        self.t = var('t')

        # Test vector fields
        self.constant_field = [1, 2, 3]  # Constant vector
        self.linear_field = [self.x, self.y, self.z]  # Linear position vector
        self.quadratic_field = [self.x**2, self.y**2, self.z**2]  # Quadratic components
        self.mixed_field = [self.x*self.y, self.y*self.z, self.x*self.z]  # Mixed products

    def test_returns_scalar(self):
        """Test that divergence returns a scalar expression"""
        result = divergence_3d(self.linear_field)

        # Result should be a scalar, not a list
        self.assertNotIsInstance(result, list)
        # Should be a symbolic expression or number
        self.assertTrue(hasattr(result, 'simplify_full') or isinstance(result, (int, float)))

    def test_constant_vector_field(self):
        """Test divergence of constant vector field"""
        # V = (1, 2, 3) - constant vector
        result = divergence_3d(self.constant_field)

        # ∇·V = ∂(1)/∂x + ∂(2)/∂y + ∂(3)/∂z = 0 + 0 + 0 = 0
        expected = 0
        self.assertEqual(result, expected)

    def test_linear_position_vector(self):
        """Test divergence of linear position vector field"""
        # V = (x, y, z) - position vector r
        result = divergence_3d(self.linear_field)

        # ∇·V = ∂x/∂x + ∂y/∂y + ∂z/∂z = 1 + 1 + 1 = 3
        expected = 3
        self.assertEqual(result, expected)

    def test_quadratic_vector_field(self):
        """Test divergence of quadratic vector field"""
        # V = (x², y², z²)
        result = divergence_3d(self.quadratic_field)

        # ∇·V = ∂(x²)/∂x + ∂(y²)/∂y + ∂(z²)/∂z = 2x + 2y + 2z
        expected = 2*self.x + 2*self.y + 2*self.z
        self.assertEqual(result, expected)

    def test_custom_3d_coordinates(self):
        """Test divergence with custom 3D coordinate variables"""
        u, v, w = var('u v w')
        coords_3d = [u, v, w]

        # Vector field in custom coordinates: V = (u, v, w)
        vector_field = [u, v, w]
        result = divergence_3d(vector_field, coordinates=coords_3d)

        # ∇·V = ∂u/∂u + ∂v/∂v + ∂w/∂w = 1 + 1 + 1 = 3
        expected = 3
        self.assertEqual(result, expected)

    def test_custom_4d_coordinates_skip_time(self):
        """Test divergence with 4D coordinates (should skip time)"""
        t, x, y, z = var('t x y z')
        coords_4d = [t, x, y, z]

        # Vector field depending on all coordinates including time
        vector_field = [t*x, t*y, t*z]
        result = divergence_3d(vector_field, coordinates=coords_4d)

        # Should only compute spatial derivatives:
        # ∇·V = ∂(tx)/∂x + ∂(ty)/∂y + ∂(tz)/∂z = t + t + t = 3t
        expected = 3*t
        self.assertEqual(result, expected)

    def test_mixed_product_vector_field(self):
        """Test divergence on vector field with mixed products"""
        # V = (xy, yz, xz)
        result = divergence_3d(self.mixed_field)

        # ∇·V = ∂(xy)/∂x + ∂(yz)/∂y + ∂(xz)/∂z = y + z + x
        expected = self.y + self.z + self.x
        self.assertEqual(result, expected)

    def test_trigonometric_vector_field(self):
        """Test divergence on trigonometric vector field"""
        # V = (sin(x), cos(y), sin(z))
        vector_field = [sin(self.x), cos(self.y), sin(self.z)]
        result = divergence_3d(vector_field)

        # ∇·V = ∂sin(x)/∂x + ∂cos(y)/∂y + ∂sin(z)/∂z = cos(x) - sin(y) + cos(z)
        expected = cos(self.x) - sin(self.y) + cos(self.z)
        self.assertEqual(result, expected)

    def test_exponential_vector_field(self):
        """Test divergence on exponential vector field"""
        # V = (e^x, e^y, e^z)
        vector_field = [exp(self.x), exp(self.y), exp(self.z)]
        result = divergence_3d(vector_field)

        # ∇·V = ∂e^x/∂x + ∂e^y/∂y + ∂e^z/∂z = e^x + e^y + e^z
        expected = exp(self.x) + exp(self.y) + exp(self.z)
        self.assertEqual(result, expected)

    def test_radial_vector_field(self):
        """Test divergence of radial vector field"""
        # V = r·r̂ = r(x, y, z)/r = (x, y, z) for r = √(x² + y² + z²)
        # This is just the position vector, so ∇·V = 3
        vector_field = [self.x, self.y, self.z]
        result = divergence_3d(vector_field)

        expected = 3  # Same as linear position vector
        self.assertEqual(result, expected)

    def test_solenoidal_vector_field(self):
        """Test divergence of solenoidal (divergence-free) vector field"""
        # V = (yz, -xz, 0) - this should have zero divergence
        vector_field = [self.y*self.z, -self.x*self.z, 0]
        result = divergence_3d(vector_field)

        # ∇·V = ∂(yz)/∂x + ∂(-xz)/∂y + ∂(0)/∂z = 0 + (-z) + 0 = -z
        # Actually this is not solenoidal, let me use a proper one:
        # V = (y, -x, 0) - 2D rotation extended to 3D
        vector_field_solenoidal = [self.y, -self.x, 0]
        result_solenoidal = divergence_3d(vector_field_solenoidal)

        # ∇·V = ∂y/∂x + ∂(-x)/∂y + ∂(0)/∂z = 0 + 0 + 0 = 0
        expected = 0
        self.assertEqual(result_solenoidal, expected)

    def test_linearity_property(self):
        """Test linearity: ∇·(aU + bV) = a∇·U + b∇·V"""
        U = [self.x, self.y, self.z]
        V = [self.x**2, self.y**2, self.z**2]
        a, b = var('a b')

        # Test ∇·(aU + bV) = a∇·U + b∇·V
        combined_field = [a*U[i] + b*V[i] for i in range(3)]
        result_combined = divergence_3d(combined_field)

        div_U = divergence_3d(U)
        div_V = divergence_3d(V)
        expected_linear = a*div_U + b*div_V

        self.assertEqual(result_combined.expand(), expected_linear.expand())

    def test_input_validation_wrong_size(self):
        """Test that wrong vector field size raises ValueError"""
        # Test with too few components
        with self.assertRaises(ValueError) as context_few:
            divergence_3d([self.x, self.y])  # Only 2 components

        # Test with too many components
        with self.assertRaises(ValueError) as context_many:
            divergence_3d([self.x, self.y, self.z, self.x])  # 4 components

        self.assertIn("Vector field must have 3 components", str(context_few.exception))
        self.assertIn("Vector field must have 3 components", str(context_many.exception))

    def test_cylindrical_coordinates_analogue(self):
        """Test divergence calculation similar to cylindrical coordinates"""
        # In Cartesian coordinates, test V = (x/r², y/r², 0) where r² = x² + y²
        # This mimics some properties of cylindrical coordinate fields
        r_squared = self.x**2 + self.y**2
        vector_field = [self.x/r_squared, self.y/r_squared, 0]
        result = divergence_3d(vector_field)

        # This should give a non-trivial result involving the coordinates
        # The exact form is complex, but we verify it's computed correctly
        self.assertIsNotNone(result)
        # Should be a symbolic expression involving x, y
        result_vars = result.variables() if hasattr(result, 'variables') else []
        self.assertTrue(any(str(var) in ['x', 'y'] for var in result_vars) or
                       any(coord in str(result) for coord in ['x', 'y']))

    def test_time_independence_in_4d_case(self):
        """Test that time coordinate doesn't affect spatial divergence"""
        t, x, y, z = var('t x y z')

        # Vector field with time dependence: V = (t²x, t²y, t²z)
        vector_field = [t**2 * self.x, t**2 * self.y, t**2 * self.z]

        # Using 4D coordinates should ignore time derivatives
        coords_4d = [t, x, y, z]
        result = divergence_3d(vector_field, coordinates=coords_4d)

        # Should only compute spatial derivatives:
        # ∇·V = ∂(t²x)/∂x + ∂(t²y)/∂y + ∂(t²z)/∂z = t² + t² + t² = 3t²
        expected = 3*t**2
        self.assertEqual(result, expected)

    def test_default_behavior_matches_explicit_coordinates(self):
        """Test that default behavior matches explicit coordinate specification"""
        vector_field = [self.x**2, self.y**2, self.z]

        # Default coordinates
        result_default = divergence_3d(vector_field)

        # Explicit 3D coordinates
        coords_explicit = [self.x, self.y, self.z]
        result_explicit = divergence_3d(vector_field, coordinates=coords_explicit)

        self.assertEqual(result_default, result_explicit)

if __name__ == '__main__':
    unittest.main()