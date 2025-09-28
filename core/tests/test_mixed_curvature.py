#!/usr/bin/env python
"""
Test suite for mixed extrinsic curvature K^i_j
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

class TestMixedExtrinsicCurvature:
    """Test mixed extrinsic curvature K^i_j"""

    def setup_method(self):
        self.M, self.X, self.t, self.r, self.th, self.ph = setup_manifold_and_coordinates()
        self.Phi, self.A, self.N = define_temporal_potential(self.M, self.t, self.r)
        self.gamma_rr, self.gamma_thth, self.gamma_phph = compute_3metric_components(self.A, self.r, self.th)
        self.K_rr, self.K_thth, self.K_phph, _, _, _ = compute_extrinsic_curvature(
            self.gamma_rr, self.gamma_thth, self.gamma_phph, self.N, self.t)

    def test_Krr_up_formula(self):
        """Test K^r_r = γ^rr K_rr = A × K_rr"""
        Krr_up, Kthth_up, Kphph_up, K_trace = compute_mixed_extrinsic_curvature(
            self.K_rr, self.K_thth, self.K_phph, self.A, self.r, self.th)

        # Test index raising formula
        expected_Krr_up = self.A * self.K_rr
        assert (Krr_up - expected_Krr_up).simplify_full() == 0, \
            "K^r_r = γ^rr K_rr = A × K_rr"

        # Verify γ^rr = A
        gamma_rr_inv = self.A
        assert (gamma_rr_inv * self.gamma_rr - 1).simplify_full() == 0, \
            "γ^rr γ_rr = 1"

        # Test dependencies
        phi_func = function('Phi')(self.t, self.r)
        assert Krr_up.has(phi_func), "K^r_r depends on Phi"
        assert Krr_up.has(self.t), "K^r_r depends on t"
        assert Krr_up.has(self.r), "K^r_r depends on r"
        assert not Krr_up.has(self.th), "K^r_r independent of θ"
        assert not Krr_up.has(self.ph), "K^r_r independent of φ"

        # Test type
        from sage.symbolic.expression import Expression
        assert isinstance(Krr_up, Expression), "K^r_r is expression"

    def test_angular_mixed_components_computed(self):
        """Test K^θ_θ and K^φ_φ computation and vanishing properties"""
        Krr_up, Kthth_up, Kphph_up, K_trace = compute_mixed_extrinsic_curvature(
            self.K_rr, self.K_thth, self.K_phph, self.A, self.r, self.th)

        # Test index raising formulas
        expected_Kthth_up = (1/self.r**2) * self.K_thth
        expected_Kphph_up = (1/(self.r**2 * sin(self.th)**2)) * self.K_phph
        assert (Kthth_up - expected_Kthth_up).simplify_full() == 0, \
            "K^θ_θ = γ^θθ K_θθ"
        assert (Kphph_up - expected_Kphph_up).simplify_full() == 0, \
            "K^φ_φ = γ^φφ K_φφ"

        # Verify preconditions and results
        assert self.K_thth == 0, "K_θθ = 0"
        assert self.K_phph == 0, "K_φφ = 0"
        assert Kthth_up == 0, "K^θ_θ = 0"
        assert Kphph_up == 0, "K^φ_φ = 0"

        # Test coordinate independence (zero should have no dependencies)
        phi_func = function('Phi')(self.t, self.r)
        assert not Kthth_up.has(phi_func), "K^θ_θ has no Phi"
        assert not Kphph_up.has(phi_func), "K^φ_φ has no Phi"
        assert not Kthth_up.has(self.t), "K^θ_θ has no t"
        assert not Kphph_up.has(self.t), "K^φ_φ has no t"

        # Test type
        from sage.symbolic.expression import Expression
        assert isinstance(Kthth_up, Expression), "K^θ_θ is expression"
        assert isinstance(Kphph_up, Expression), "K^φ_φ is expression"

    def test_trace_formula(self):
        """Test K = K^i_i = K^r_r and critical trace simplification for flux law"""
        Krr_up, Kthth_up, Kphph_up, K_trace = compute_mixed_extrinsic_curvature(
            self.K_rr, self.K_thth, self.K_phph, self.A, self.r, self.th)

        # Test trace formula
        expected_trace = Krr_up + Kthth_up + Kphph_up
        assert (K_trace - expected_trace).simplify_full() == 0, \
            "K = K^r_r + K^θ_θ + K^φ_φ"

        # Verify angular components vanish
        assert Kthth_up == 0, "K^θ_θ = 0"
        assert Kphph_up == 0, "K^φ_φ = 0"

        # Test critical simplification: K = K^r_r
        assert (K_trace - Krr_up).simplify_full() == 0, "K = K^r_r"

        # Test coordinate dependencies
        phi_func = function('Phi')(self.t, self.r)
        assert K_trace.has(phi_func), "K depends on Phi"
        assert K_trace.has(self.t), "K depends on t"
        assert K_trace.has(self.r), "K depends on r"
        assert not K_trace.has(self.th), "K independent of θ"
        assert not K_trace.has(self.ph), "K independent of φ"

        # Verify Y^r_r = K^r_r - K = 0
        Y_rr = Krr_up - K_trace
        assert Y_rr.simplify_full() == 0, "Y^r_r = 0"

        # Verify Y^θ_θ = -K and Y^φ_φ = -K
        Y_thth = Kthth_up - K_trace
        Y_phph = Kphph_up - K_trace
        assert (Y_thth + K_trace).simplify_full() == 0, "Y^θ_θ = -K"
        assert (Y_phph + K_trace).simplify_full() == 0, "Y^φ_φ = -K"

        # Test type
        from sage.symbolic.expression import Expression
        assert isinstance(K_trace, Expression), f"K is symbolic expression"

    def test_Krr_up_simplified_form(self):
        """Test K^r_r simplifies to e^(-Φ) ∂_tΦ - key formula for flux law"""
        Krr_up, Kthth_up, Kphph_up, K_trace = compute_mixed_extrinsic_curvature(
            self.K_rr, self.K_thth, self.K_phph, self.A, self.r, self.th)

        # Test the critical simplification
        Krr_up_simplified = Krr_up.simplify_full()
        expected = exp(-self.Phi.expr()) * diff(self.Phi.expr(), self.t)
        assert (Krr_up_simplified - expected).simplify_full() == 0, \
            "K^r_r should simplify to exp(-Φ) ∂_t Φ"

        # Verify the mathematical derivation step by step
        # K_rr = -(1/2N) ∂_t γ_rr
        dt_gamma_rr = diff(self.gamma_rr, self.t)
        K_rr_check = -(1/(2*self.N)) * dt_gamma_rr
        assert (self.K_rr - K_rr_check).simplify_full() == 0, "K_rr formula check"

        # K^r_r = A × K_rr
        Krr_up_check = self.A * self.K_rr
        assert (Krr_up - Krr_up_check).simplify_full() == 0, "K^r_r = A × K_rr"

        # Verify γ_rr = 1/A = exp(-2*Phi)
        gamma_rr_check = exp(-2*self.Phi.expr())
        assert (self.gamma_rr - gamma_rr_check).simplify_full() == 0, "γ_rr = exp(-2Φ)"

        # Verify ∂_t γ_rr = -2*exp(-2*Phi)*∂_t Phi
        dt_gamma_rr_check = -2*exp(-2*self.Phi.expr())*diff(self.Phi.expr(), self.t)
        assert (dt_gamma_rr - dt_gamma_rr_check).simplify_full() == 0, \
            "∂_t γ_rr = -2*exp(-2Φ)*∂_t Φ"

        # Verify trace equals K^r_r
        assert (K_trace - Krr_up).simplify_full() == 0, "K = K^r_r"
        assert (K_trace - Krr_up_simplified).simplify_full() == 0, \
            "K = exp(-Φ) ∂_t Φ"

        # Verify Y^r_r = 0 using simplified form
        Y_rr = Krr_up_simplified - K_trace
        assert Y_rr.simplify_full() == 0, "Y^r_r = K^r_r - K = 0"

        # Test coordinate dependencies
        phi_func = function('Phi')(self.t, self.r)
        assert Krr_up_simplified.has(phi_func), "K^r_r contains Phi"
        assert Krr_up_simplified.has(self.t), "K^r_r depends on t"
        assert Krr_up_simplified.has(self.r), "K^r_r depends on r"
        assert not Krr_up_simplified.has(self.th), "K^r_r independent of θ"
        assert not Krr_up_simplified.has(self.ph), "K^r_r independent of φ"

        # Verify the form can be written as (1/N) ∂_t Phi
        alternative_form = (1/self.N) * diff(self.Phi.expr(), self.t)
        assert (Krr_up_simplified - alternative_form).simplify_full() == 0, \
            "K^r_r = (1/N) ∂_t Φ = exp(-Φ) ∂_t Φ"

if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
