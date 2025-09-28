#!/usr/bin/env python
"""
Test suite for Y tensor Y^j_i = K^j_i - δ^j_i K
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

class TestYTensor:
    """Test Y tensor Y^j_i = K^j_i - δ^j_i K"""

    def test_Y_definition_verification(self):
        """Test Y^j_i = K^j_i - δ^j_i K definition"""
        Y_rr, Y_thth, Y_phph = compute_Y_tensor(self.Krr_up, self.K_trace)

        assert (Y_rr - (self.Krr_up - self.K_trace)).simplify_full() == 0
        assert Y_thth == -self.K_trace
        assert Y_phph == -self.K_trace

        assert (self.Kthth_up).simplify_full() == 0
        assert (self.Kphph_up).simplify_full() == 0

    def setup_method(self):
        self.M, self.X, self.t, self.r, self.th, self.ph = setup_manifold_and_coordinates()
        self.Phi, self.A, self.N = define_temporal_potential(self.M, self.t, self.r)
        self.gamma_rr, self.gamma_thth, self.gamma_phph = compute_3metric_components(self.A, self.r, self.th)
        self.K_rr, self.K_thth, self.K_phph, _, _, _ = compute_extrinsic_curvature(
            self.gamma_rr, self.gamma_thth, self.gamma_phph, self.N, self.t)
        self.Krr_up, self.Kthth_up, self.Kphph_up, self.K_trace = compute_mixed_extrinsic_curvature(
            self.K_rr, self.K_thth, self.K_phph, self.A, self.r, self.th)

    def test_Y_rr_zero(self):
        """Test Y^r_r = K^r_r - K = 0 (crucial property for flux law derivation)"""
        Y_rr, Y_thth, Y_phph = compute_Y_tensor(self.Krr_up, self.K_trace)

        expected_Y_rr = self.Krr_up - self.K_trace
        assert (Y_rr - expected_Y_rr).simplify_full() == 0

        assert Y_rr.simplify_full() == 0
        assert (self.K_trace - self.Krr_up).simplify_full() == 0

    def test_Y_angular_components(self):
        """Test Y^θ_θ = Y^φ_φ = -K"""
        Y_rr, Y_thth, Y_phph = compute_Y_tensor(self.Krr_up, self.K_trace)

        assert Y_thth == -self.K_trace
        assert Y_phph == -self.K_trace
        assert (Y_thth - Y_phph).simplify_full() == 0

    def test_Y_tensor_trace_calculation(self):
        """Test Y tensor trace calculation"""
        Y_rr, Y_thth, Y_phph = compute_Y_tensor(self.Krr_up, self.K_trace)

        Y_trace = Y_rr + Y_thth + Y_phph
        expected_trace = 0 + (-self.K_trace) + (-self.K_trace)
        assert (Y_trace - expected_trace).simplify_full() == 0

        substituted_trace = -2*self.K_trace
        assert (Y_trace - substituted_trace).simplify_full() == 0

if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
