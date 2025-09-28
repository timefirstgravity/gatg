#!/usr/bin/env python
"""
Test suite for extrinsic curvature computation
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
    from core import compute_extrinsic_curvature, compute_mixed_extrinsic_curvature, compute_Y_tensor, compute_covariant_divergence
except ImportError as e:
    pytest.skip(f"Sage not available: {e}", allow_module_level=True)

class TestExtrinsicCurvature:
    """Test extrinsic curvature computation"""

    def setup_method(self):
        self.M, self.X, self.t, self.r, self.th, self.ph = setup_manifold_and_coordinates()
        self.Phi, self.A, self.N = define_temporal_potential(self.M, self.t, self.r)
        self.gamma_rr, self.gamma_thth, self.gamma_phph = compute_3metric_components(self.A, self.r, self.th)

    def test_time_derivatives(self):
        """Test time derivatives of 3-metric components"""
        K_rr, K_thth, K_phph, dt_gamma_rr, dt_gamma_thth, dt_gamma_phph = compute_extrinsic_curvature(
            self.gamma_rr, self.gamma_thth, self.gamma_phph, self.N, self.t)

        # Test types
        from sage.symbolic.expression import Expression
        assert isinstance(dt_gamma_rr, Expression), "∂_t γ_rr is expression"
        assert isinstance(dt_gamma_thth, Expression), "∂_t γ_θθ is expression"
        assert isinstance(dt_gamma_phph, Expression), "∂_t γ_φφ is expression"

        # Test correctness
        expected_dt_gamma_rr = diff(1/self.A, self.t)
        assert (dt_gamma_rr - expected_dt_gamma_rr).simplify_full() == 0, "∂_t γ_rr = ∂_t(1/A)"

        # Test angular components vanish
        assert dt_gamma_thth == 0, "∂_t γ_θθ = 0"
        assert dt_gamma_phph == 0, "∂_t γ_φφ = 0"

        # Test dependencies
        phi_func = function('Phi')(self.t, self.r)
        assert dt_gamma_rr.has(phi_func), "∂_t γ_rr contains Phi"
        assert dt_gamma_rr.has(self.t), "∂_t γ_rr depends on t"
        assert dt_gamma_rr.has(self.r), "∂_t γ_rr depends on r"
        assert not dt_gamma_rr.has(self.th), "∂_t γ_rr independent of θ"
        assert not dt_gamma_rr.has(self.ph), "∂_t γ_rr independent of φ"

        # Check derivative structure exists
        str_repr = str(dt_gamma_rr)
        assert 'diff' in str_repr or 'D[' in str_repr, "∂_t γ_rr has derivative"

    def test_K_rr_formula(self):
        """Test K_rr = -(1/2N) ∂_t γ_rr"""
        K_rr, K_thth, K_phph, dt_gamma_rr, dt_gamma_thth, dt_gamma_phph = compute_extrinsic_curvature(
            self.gamma_rr, self.gamma_thth, self.gamma_phph, self.N, self.t)

        # Test ADM formula
        expected_K_rr = -(1/(2*self.N)) * dt_gamma_rr
        assert (K_rr - expected_K_rr).simplify_full() == 0, "K_rr = -(1/2N) ∂_t γ_rr"

        # Test dependencies
        phi_func = function('Phi')(self.t, self.r)
        assert K_rr.has(phi_func), "K_rr depends on Phi"
        assert K_rr.has(self.t), "K_rr depends on t"
        assert K_rr.has(self.r), "K_rr depends on r"

        # Test type
        from sage.symbolic.expression import Expression
        assert isinstance(K_rr, Expression), "K_rr is expression"

    def test_angular_K_components_zero(self):
        """Test K_θθ = K_φφ = 0"""
        K_rr, K_thth, K_phph, dt_gamma_rr, dt_gamma_thth, dt_gamma_phph = compute_extrinsic_curvature(
            self.gamma_rr, self.gamma_thth, self.gamma_phph, self.N, self.t)

        # Test vanishing components
        assert K_thth == 0, "K_θθ = 0"
        assert K_phph == 0, "K_φφ = 0"

        # Test formula consistency
        expected_K_thth = -(1/(2*self.N)) * dt_gamma_thth
        expected_K_phph = -(1/(2*self.N)) * dt_gamma_phph
        assert (K_thth - expected_K_thth).simplify_full() == 0, "K_θθ = -(1/2N) ∂_t γ_θθ"
        assert (K_phph - expected_K_phph).simplify_full() == 0, "K_φφ = -(1/2N) ∂_t γ_φφ"

        # Verify time derivatives vanish
        assert dt_gamma_thth == 0, "∂_t γ_θθ = 0"
        assert dt_gamma_phph == 0, "∂_t γ_φφ = 0"

        # Test coordinate independence
        phi_func = function('Phi')(self.t, self.r)
        assert not K_thth.has(phi_func), "K_θθ has no Phi"
        assert not K_phph.has(phi_func), "K_φφ has no Phi"
        assert not K_thth.has(self.t), "K_θθ has no t"
        assert not K_phph.has(self.t), "K_φφ has no t"

        # Test type
        from sage.symbolic.expression import Expression
        assert isinstance(K_thth, Expression), "K_θθ is expression"
        assert isinstance(K_phph, Expression), "K_φφ is expression"

if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
