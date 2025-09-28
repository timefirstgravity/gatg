#!/usr/bin/env python
"""
Test suite for covariant divergence computation
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

class TestCovariantDivergence:
    """Test covariant divergence computation"""

    def test_Y_component_symmetry(self):
        """Test Y component symmetry in spherical coordinates"""
        DY, Gamma_theta_rtheta, Gamma_phi_rphi = compute_covariant_divergence(
            self.Y_thth, self.Y_phph, self.r)

        assert (self.Y_thth - self.Y_phph).simplify_full() == 0

        expected_DY = -2 * Gamma_theta_rtheta * self.Y_thth
        assert (DY - expected_DY).simplify_full() == 0

    def setup_method(self):
        # Setup core components only
        M, X, t, r, th, ph = setup_manifold_and_coordinates()
        Phi, A, N = define_temporal_potential(M, t, r)
        gamma_rr, gamma_thth, gamma_phph = compute_3metric_components(A, r, th)
        K_rr, K_thth, K_phph, dt_gamma_rr, dt_gamma_thth, dt_gamma_phph = compute_extrinsic_curvature(
            gamma_rr, gamma_thth, gamma_phph, N, t)
        Krr_up, Kthth_up, Kphph_up, K_trace = compute_mixed_extrinsic_curvature(
            K_rr, K_thth, K_phph, A, r, th)
        self.Y_rr, self.Y_thth, self.Y_phph = compute_Y_tensor(Krr_up, K_trace)
        self.r = r

    def test_christoffel_symbols(self):
        """Test Christoffel symbols for spherical metric"""
        DY, Gamma_theta_rtheta, Gamma_phi_rphi = compute_covariant_divergence(
            self.Y_thth, self.Y_phph, self.r)

        assert (Gamma_theta_rtheta - 1/self.r).simplify_full() == 0
        assert (Gamma_phi_rphi - 1/self.r).simplify_full() == 0

        assert Gamma_theta_rtheta == Gamma_phi_rphi

        expected_gamma = 1/self.r
        assert (Gamma_theta_rtheta - expected_gamma).simplify_full() == 0

    def test_divergence_formula(self):
        """Test D_j Y^j_r = -Γ^θ_rθ Y^θ_θ - Γ^φ_rφ Y^φ_φ"""
        DY, Gamma_theta_rtheta, Gamma_phi_rphi = compute_covariant_divergence(
            self.Y_thth, self.Y_phph, self.r)

        expected_DY = -Gamma_theta_rtheta * self.Y_thth - Gamma_phi_rphi * self.Y_phph
        assert (DY - expected_DY).simplify_full() == 0

        since_Y_components_equal = (self.Y_thth - self.Y_phph).simplify_full() == 0
        if since_Y_components_equal:
            simplified_expected = -2 * Gamma_theta_rtheta * self.Y_thth
            assert (DY - simplified_expected).simplify_full() == 0

    def test_divergence_simplified_form(self):
        """Test D_j Y^j_r mathematical structure"""
        DY, Gamma_theta_rtheta, Gamma_phi_rphi = compute_covariant_divergence(
            self.Y_thth, self.Y_phph, self.r)

        DY_simplified = DY.simplify_full()

        manual_calculation = -2 * (1/self.r) * self.Y_thth
        assert (DY_simplified - manual_calculation).simplify_full() == 0

        assert DY_simplified.has(1/self.r)

if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
