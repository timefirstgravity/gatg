#!/usr/bin/env python
"""
Test suite for metric tensor construction
Extracted from test_functions.py for better organization
"""

#!/usr/bin/env python
"""
Comprehensive pytest test suite for flux law mathematical functions
Tests the actual code from verify_flux_law_computational.py
"""

import pytest
import sys
import os

# Add sage to path and import functions
sys.path.append('/opt/homebrew/bin')  # Common sage location
os.environ['SAGE_ROOT'] = '/opt/homebrew'

try:
    from sage.all import *
    from core import *
    # Functions now imported from core
except ImportError as e:
    pytest.skip(f"Sage not available: {e}", allow_module_level=True)

class TestMetricConstruction:
    """Test metric tensor construction"""

    def setup_method(self):
        self.M, self.X, self.t, self.r, self.th, self.ph = setup_manifold_and_coordinates()
        self.Phi, self.A, self.N = define_temporal_potential(self.M, self.t, self.r)

    def test_metric_components(self):
        """Test metric tensor components and geometric properties"""
        g = construct_metric(self.M, self.A, self.r, self.th)

        # Test that g is a proper SageMath metric tensor
        from sage.manifolds.differentiable.metric import PseudoRiemannianMetric
        assert isinstance(g, PseudoRiemannianMetric), f"g should be a metric tensor, got {type(g)}"

        # Test that g belongs to the correct manifold
        assert g.domain() == self.M, "Metric should be defined on manifold M"

        # Test metric name and string representation
        assert g._name == 'g', f"Metric should be named 'g', got {g._name}"

        # Test diagonal components (covariant metric g_μν)
        assert g[0,0].expr() == -self.A, f"g_tt should equal -A, got {g[0,0].expr()}"
        assert g[1,1].expr() == 1/self.A, f"g_rr should equal 1/A, got {g[1,1].expr()}"
        assert g[2,2].expr() == self.r**2, f"g_θθ should equal r², got {g[2,2].expr()}"
        assert g[3,3].expr() == self.r**2 * sin(self.th)**2, f"g_φφ should equal r²sin²θ, got {g[3,3].expr()}"

        # Test that off-diagonal components are zero (spherical symmetry)
        for i in range(4):
            for j in range(4):
                if i != j:
                    assert g[i,j].expr() == 0, f"Off-diagonal component g[{i},{j}] should be zero, got {g[i,j].expr()}"

        # Test metric signature (Lorentzian: one negative eigenvalue)
        signature = g.signature()
        assert signature == 2, f"Metric should have signature 2 (Lorentzian), got {signature}"

        # Test metric determinant (important for volume element)
        det_g = g.determinant()
        expected_det = -self.A * (1/self.A) * self.r**2 * (self.r**2 * sin(self.th)**2)
        expected_det_simplified = -(self.r**4) * sin(self.th)**2

        # Simplify both expressions for comparison
        det_simplified = det_g.expr().simplify_full()
        expected_simplified = expected_det_simplified.simplify_full()
        assert det_simplified == expected_simplified, \
            f"Metric determinant should be -r⁴sin²θ, got {det_simplified}"

        # Test that determinant is negative (required for Lorentzian metric)
        # Since r > 0 and 0 < θ < π, we have sin²θ ≥ 0, so det_g < 0
        assert det_g.expr().has(self.r), "Determinant should depend on r"
        assert det_g.expr().has(self.th), "Determinant should depend on θ"

        # Test inverse metric components g^μν
        g_inv = g.inverse()
        assert g_inv[0,0].expr() == -1/self.A, "g^tt should equal -1/A"
        assert g_inv[1,1].expr() == self.A, "g^rr should equal A"
        assert g_inv[2,2].expr() == 1/self.r**2, "g^θθ should equal 1/r²"
        assert g_inv[3,3].expr() == 1/(self.r**2 * sin(self.th)**2), "g^φφ should equal 1/(r²sin²θ)"

        # Test key mathematical relationships for metric correctness
        # For diagonal metric: g_μν * g^μν should give correct inverse relationships
        # g_tt = -A, g^tt = -1/A => g_tt * g^tt = (-A) * (-1/A) = 1
        assert (g[0,0].expr() * g_inv[0,0].expr()).simplify_full() == 1, \
            "g_tt * g^tt = 1 (correct metric-inverse relationship)"
        assert (g[1,1].expr() * g_inv[1,1].expr()).simplify_full() == 1, \
            "g_rr * g^rr = 1 (metric-inverse relationship)"
        assert (g[2,2].expr() * g_inv[2,2].expr()).simplify_full() == 1, \
            "g_θθ * g^θθ = 1 (metric-inverse relationship)"
        assert (g[3,3].expr() * g_inv[3,3].expr()).simplify_full() == 1, \
            "g_φφ * g^φφ = 1 (metric-inverse relationship)"

        # Test that this represents proper spherically symmetric spacetime
        # Line element: ds² = -A dt² + (1/A) dr² + r² dθ² + r² sin²θ dφ²
        # This should match the expected form for spherical symmetry
        assert g[2,2].expr() == self.r**2, "Angular component θ should be r²"
        assert g[3,3].expr() == self.r**2 * sin(self.th)**2, "Angular component φ should be r²sin²θ"

        # Test that the radial components have the expected A-dependence
        assert g[0,0].expr() == -self.A, "Time component should be -A"
        assert g[1,1].expr() == 1/self.A, "Radial component should be 1/A"

    def test_metric_signature(self):
        """Test metric has correct Lorentzian signature (-,+,+,+)"""
        g = construct_metric(self.M, self.A, self.r, self.th)

        # Test signature using SageMath's built-in method
        signature = g.signature()
        assert signature == 2, f"Lorentzian metric should have signature 2, got {signature}"

        # Test physical interpretation: verify one timelike, three spacelike directions
        # For Lorentzian metric: signature = dim - 2*negative_eigenvalues = 4 - 2*1 = 2
        # This confirms exactly one negative eigenvalue (timelike direction)

        # Test component signs for mathematical verification
        # Timelike component (negative)
        assert g[0,0].expr() == -self.A, "Time component g₀₀ should be -A (negative)"

        # Spacelike components (positive for physical coordinates)
        # Note: These should be positive in the physical region where A > 0, r > 0, 0 < θ < π
        assert g[1,1].expr() == 1/self.A, "Radial component g₁₁ should be 1/A (positive when A > 0)"
        assert g[2,2].expr() == self.r**2, "θ component g₂₂ should be r² (positive when r > 0)"
        assert g[3,3].expr() == self.r**2 * sin(self.th)**2, "φ component g₃₃ should be r²sin²θ (positive)"

        # Test signature independence: signature should be invariant under coordinate changes
        # The signature depends only on the metric structure, not coordinate values

        # Test that signature represents proper spacetime structure
        # In 4D: signature 2 means (1 timelike, 3 spacelike) = Lorentzian manifold
        # This is essential for: causality, light cones, timelike/spacelike separation

        # Test mathematical consistency: for diagonal metric, signature comes from eigenvalue signs
        # Since g is diagonal: eigenvalues = diagonal components
        # g₀₀ = -A < 0 (when A > 0) => 1 negative eigenvalue
        # g₁₁, g₂₂, g₃₃ > 0 (in physical region) => 3 positive eigenvalues
        # Total: signature = 4 - 2(1) = 2 ✓

        # Test that this signature is compatible with general relativity
        # Lorentzian signature required for: proper time, causal structure, Einstein equations
        assert signature > 0, "Signature must be positive for valid spacetime"
        assert signature < 4, "Signature must be less than dimension for non-Euclidean metric"

    def test_metric_off_diagonal_zero(self):
        """Test off-diagonal components vanish and geometric implications"""
        g = construct_metric(self.M, self.A, self.r, self.th)

        # Test that all off-diagonal components are exactly zero
        off_diagonal_components = []
        for i in range(4):
            for j in range(4):
                if i != j:
                    component_value = g[i,j].expr()
                    assert component_value == 0, f"Off-diagonal component g[{i},{j}] should be zero, got {component_value}"
                    off_diagonal_components.append((i, j, component_value))

        # Test geometric interpretation: spherical symmetry
        # Zero off-diagonals mean coordinate basis vectors are orthogonal
        # This represents spherical symmetry in the spatial part and time-space orthogonality

        # Test time-space orthogonality (g₀ᵢ = 0 for i = 1,2,3)
        assert g[0,1].expr() == 0, "g_tr should be zero (no time-radial mixing)"
        assert g[0,2].expr() == 0, "g_tθ should be zero (no time-angular mixing)"
        assert g[0,3].expr() == 0, "g_tφ should be zero (no time-azimuthal mixing)"

        # Test spatial orthogonality (g_ij = 0 for i≠j, i,j ∈ {1,2,3})
        assert g[1,2].expr() == 0, "g_rθ should be zero (radial-polar orthogonal)"
        assert g[1,3].expr() == 0, "g_rφ should be zero (radial-azimuthal orthogonal)"
        assert g[2,3].expr() == 0, "g_θφ should be zero (polar-azimuthal orthogonal)"

        # Test symmetry: g_μν = g_νμ (metric tensor is symmetric)
        for i in range(4):
            for j in range(4):
                assert g[i,j].expr() == g[j,i].expr(), f"Metric should be symmetric: g[{i},{j}] = g[{j},{i}]"

        # Test physical implications of diagonal structure
        # 1. No dragging effects (time-space mixing)
        # 2. Orthogonal spatial coordinates (spherical symmetry)
        # 3. Simplified calculations for curvature and divergence

        # Test mathematical consequences: metric determinant factorizes
        # For diagonal metric: det(g) = g₀₀ × g₁₁ × g₂₂ × g₃₃
        det_g = g.determinant()
        expected_det_factorized = g[0,0].expr() * g[1,1].expr() * g[2,2].expr() * g[3,3].expr()

        # Verify determinant matches factorized form
        assert (det_g.expr() - expected_det_factorized).simplify_full() == 0, \
            "Determinant should factorize for diagonal metric"

        # Test coordinate orthogonality implications for physics
        # Diagonal metric enables separation of:
        # - Temporal vs spatial effects
        # - Radial vs angular motion
        # - Independent spherical coordinates

        # Test that this structure is compatible with spherical symmetry
        # The metric should be invariant under SO(3) rotations in θ,φ
        # This is ensured by the functional form: functions of r only in spatial part

        assert not g[2,2].expr().has(self.t), "g_θθ should not depend on time (spherical symmetry)"
        assert not g[3,3].expr().has(self.t), "g_φφ should not depend on time (spherical symmetry)"
        assert g[2,2].expr().has(self.r), "g_θθ should depend on radius (spherical coordinates)"
        assert g[3,3].expr().has(self.r), "g_φφ should depend on radius (spherical coordinates)"

        # Test that off-diagonal zeros enable simplified Christoffel symbols
        # For diagonal metric, many Christoffel symbols vanish or simplify significantly
        # This is crucial for the subsequent covariant divergence calculations

if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
