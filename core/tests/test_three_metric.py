#!/usr/bin/env python
"""
Test suite for 3-metric components
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

class TestThreeMetric:
    """Test 3-metric components"""

    def setup_method(self):
        self.M, self.X, self.t, self.r, self.th, self.ph = setup_manifold_and_coordinates()
        self.Phi, self.A, self.N = define_temporal_potential(self.M, self.t, self.r)

    def test_3metric_components(self):
        """Test spatial 3-metric components"""
        gamma_rr, gamma_thth, gamma_phph = compute_3metric_components(self.A, self.r, self.th)

        # Test component values
        assert gamma_rr == 1/self.A, "γ_rr = 1/A"
        assert gamma_thth == self.r**2, "γ_θθ = r²"
        assert gamma_phph == self.r**2 * sin(self.th)**2, "γ_φφ = r²sin²θ"

        # Test types
        from sage.symbolic.expression import Expression
        assert isinstance(gamma_rr, Expression), "γ_rr is expression"
        assert isinstance(gamma_thth, Expression), "γ_θθ is expression"
        assert isinstance(gamma_phph, Expression), "γ_φφ is expression"

        # Test dependencies for γ_rr
        phi_func = function('Phi')(self.t, self.r)
        assert gamma_rr.has(phi_func), "γ_rr depends on Phi"
        assert gamma_rr.has(self.t), "γ_rr depends on t"
        assert gamma_rr.has(self.r), "γ_rr depends on r"
        assert not gamma_rr.has(self.th), "γ_rr independent of θ"

        # Test dependencies for γ_θθ
        assert gamma_thth.has(self.r), "γ_θθ depends on r"
        assert not gamma_thth.has(phi_func), "γ_θθ independent of Phi"
        assert not gamma_thth.has(self.t), "γ_θθ independent of t"
        assert not gamma_thth.has(self.th), "γ_θθ independent of θ"

        # Test dependencies for γ_φφ
        assert gamma_phph.has(self.r), "γ_φφ depends on r"
        assert gamma_phph.has(self.th), "γ_φφ depends on θ"
        assert not gamma_phph.has(phi_func), "γ_φφ independent of Phi"
        assert not gamma_phph.has(self.t), "γ_φφ independent of t"

        # Test 3-metric determinant
        det_gamma = gamma_rr * gamma_thth * gamma_phph
        expected_det = (self.r**4 * sin(self.th)**2) / self.A
        assert (det_gamma - expected_det).simplify_full() == 0, "det(γ) = r⁴sin²θ/A"

        # Test inverse relationships
        gamma_rr_inv = self.A
        gamma_thth_inv = 1/self.r**2
        gamma_phph_inv = 1/(self.r**2 * sin(self.th)**2)

        assert (gamma_rr * gamma_rr_inv).simplify_full() == 1, "γ_rr × γ^rr = 1"
        assert (gamma_thth * gamma_thth_inv).simplify_full() == 1, "γ_θθ × γ^θθ = 1"
        assert (gamma_phph * gamma_phph_inv).simplify_full() == 1, "γ_φφ × γ^φφ = 1"

    def test_3metric_inverse_relation(self):
        """Test 3-metric inverse tensor properties"""
        gamma_rr, gamma_thth, gamma_phph = compute_3metric_components(self.A, self.r, self.th)

        # Define inverse components
        gamma_rr_inv = self.A
        gamma_thth_inv = 1/self.r**2
        gamma_phph_inv = 1/(self.r**2 * sin(self.th)**2)

        # Test inverse relationships
        assert (gamma_rr * gamma_rr_inv).simplify_full() == 1, "γ_rr × γ^rr = 1"
        assert (gamma_thth * gamma_thth_inv).simplify_full() == 1, "γ_θθ × γ^θθ = 1"
        assert (gamma_phph * gamma_phph_inv).simplify_full() == 1, "γ_φφ × γ^φφ = 1"

        # Test determinant properties
        det_gamma = gamma_rr * gamma_thth * gamma_phph
        det_gamma_inv = gamma_rr_inv * gamma_thth_inv * gamma_phph_inv
        assert (det_gamma * det_gamma_inv).simplify_full() == 1, "det(γ) × det(γ⁻¹) = 1"
        assert (det_gamma_inv - 1/det_gamma).simplify_full() == 0, "det(γ⁻¹) = 1/det(γ)"

        # Test index raising and lowering
        V_r = 1
        V_th = self.r**2
        V_ph = self.r**2 * sin(self.th)**2

        V_r_raised = gamma_rr_inv * V_r
        V_th_raised = gamma_thth_inv * V_th
        V_ph_raised = gamma_phph_inv * V_ph

        V_r_lowered = gamma_rr * V_r_raised
        V_th_lowered = gamma_thth * V_th_raised
        V_ph_lowered = gamma_phph * V_ph_raised

        assert (V_r_lowered - V_r).simplify_full() == 0, "Raise then lower r"
        assert (V_th_lowered - V_th).simplify_full() == 0, "Raise then lower θ"
        assert (V_ph_lowered - V_ph).simplify_full() == 0, "Raise then lower φ"

        # Test dependencies of inverse components
        phi_func = function('Phi')(self.t, self.r)
        assert gamma_rr_inv.has(phi_func), "γ^rr depends on Phi"
        assert gamma_thth_inv.has(self.r), "γ^θθ depends on r"
        assert gamma_phph_inv.has(self.r), "γ^φφ depends on r"
        assert gamma_phph_inv.has(self.th), "γ^φφ depends on θ"

        # Test double inverse
        assert (1/gamma_rr_inv - gamma_rr).simplify_full() == 0, "(γ^rr)⁻¹ = γ_rr"
        assert (1/gamma_thth_inv - gamma_thth).simplify_full() == 0, "(γ^θθ)⁻¹ = γ_θθ"
        assert (1/gamma_phph_inv - gamma_phph).simplify_full() == 0, "(γ^φφ)⁻¹ = γ_φφ"

if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
