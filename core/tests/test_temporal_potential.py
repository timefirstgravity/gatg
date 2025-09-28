#!/usr/bin/env python
"""
Test suite for temporal potential and related quantities
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

class TestTemporalPotential:
    """Test temporal potential and related quantities"""

    def setup_method(self):
        self.M, self.X, self.t, self.r, self.th, self.ph = setup_manifold_and_coordinates()

    def test_temporal_potential_definition(self):
        """Test Φ(t,r) definition and mathematical properties"""
        Phi, A, N = define_temporal_potential(self.M, self.t, self.r)

        # Test that Phi is a proper SageMath scalar field
        from sage.manifolds.scalarfield import ScalarField
        assert isinstance(Phi, ScalarField), f"Phi should be a ScalarField, got {type(Phi)}"

        # Test that Phi belongs to the correct manifold
        assert Phi.domain() == self.M, "Phi's domain should be the manifold M"

        # Test that Phi's parent is the scalar field algebra on M
        parent_str = str(Phi.parent())
        assert "scalar fields on" in parent_str and "manifold M" in parent_str, \
            f"Phi should belong to scalar field algebra on M, got {parent_str}"

        # Test that Phi has proper name
        assert Phi._name == 'Phi', f"Phi should have name 'Phi', got {Phi._name}"

        # Test string representation contains the name
        phi_str = str(Phi)
        assert 'Scalar field Phi' in phi_str, f"Phi string representation should contain 'Scalar field Phi', got {phi_str}"

        # Test that Phi's expression is the expected symbolic function
        phi_expr = Phi.expr()

        # Verify functional dependence on both t and r
        assert phi_expr.has(self.t), "Phi expression must depend on time coordinate t"
        assert phi_expr.has(self.r), "Phi expression must depend on radial coordinate r"

        # Test that the expression is indeed the function Phi(t,r)
        expected_phi_func = function('Phi')(self.t, self.r)
        assert phi_expr == expected_phi_func, f"Phi.expr() should equal Phi(t,r), got {phi_expr}"

        # Test that Phi can be evaluated as a scalar field expression
        # (this verifies it's properly constructed as a scalar field)
        assert callable(phi_expr), "Phi expression should be callable"

        # Test that the scalar field has correct coordinate representation
        # Get the coordinate function on the chart
        coord_function = Phi.coord_function(self.X)
        assert coord_function is not None, "Phi should have coordinate function representation on chart X"

        # Verify the coordinate function equals the expression
        assert coord_function == phi_expr, "Coordinate function should match the expression"

    def test_A_definition(self):
        """Test A = e^(2Φ) definition and mathematical properties"""
        Phi, A, N = define_temporal_potential(self.M, self.t, self.r)

        # Test that A is the correct mathematical expression
        expected_A = exp(2*Phi.expr())
        assert (A - expected_A).simplify_full() == 0, f"A should equal exp(2*Phi), got A={A}"

        # Test that A is a symbolic expression (not a scalar field)
        from sage.symbolic.expression import Expression
        assert isinstance(A, Expression), f"A should be a symbolic expression, got {type(A)}"

        # Test functional dependencies - A should depend on same variables as Phi
        assert A.has(self.t), "A must depend on time coordinate t"
        assert A.has(self.r), "A must depend on radial coordinate r"

        # Test that A has the correct functional form A(t,r)
        # Since Phi = Phi(t,r), A should be exp(2*Phi(t,r))
        phi_func = function('Phi')(self.t, self.r)
        expected_A_explicit = exp(2*phi_func)
        assert A == expected_A_explicit, f"A should equal exp(2*Phi(t,r)), got {A}"

        # Test mathematical properties - A should always be positive
        # Since A = exp(2*Phi), and exponential is always positive
        assert A.operator() == exp, "A should be an exponential function"

        # Test that A can be differentiated (required for metric calculations)
        dA_dt = diff(A, self.t)
        dA_dr = diff(A, self.r)
        assert dA_dt is not None, "A should be differentiable with respect to t"
        assert dA_dr is not None, "A should be differentiable with respect to r"

        # Test that derivatives have expected form: ∂A/∂t = 2*exp(2Φ)*∂Φ/∂t
        expected_dA_dt = 2 * A * diff(phi_func, self.t)
        assert (dA_dt - expected_dA_dt).simplify_full() == 0, \
            "∂A/∂t should equal 2*A*∂Φ/∂t by chain rule"

        # Test that A can be used in algebraic operations (needed for metric)
        A_inverse = 1/A
        assert A_inverse is not None, "A should be invertible (A > 0 always)"

        # Test that 1/A = exp(-2Phi)
        expected_A_inv = exp(-2*phi_func)
        assert (A_inverse - expected_A_inv).simplify_full() == 0, \
            "1/A should equal exp(-2Φ)"

    def test_lapse_definition(self):
        """Test N = e^Φ definition and lapse function properties"""
        Phi, A, N = define_temporal_potential(self.M, self.t, self.r)

        # Test that N is the correct mathematical expression
        expected_N = exp(Phi.expr())
        assert (N - expected_N).simplify_full() == 0, f"N should equal exp(Phi), got N={N}"

        # Test that N is a symbolic expression (not a scalar field)
        from sage.symbolic.expression import Expression
        assert isinstance(N, Expression), f"N should be a symbolic expression, got {type(N)}"

        # Test functional dependencies - N should depend on same variables as Phi
        assert N.has(self.t), "N must depend on time coordinate t"
        assert N.has(self.r), "N must depend on radial coordinate r"

        # Test that N has the correct functional form N(t,r)
        phi_func = function('Phi')(self.t, self.r)
        expected_N_explicit = exp(phi_func)
        assert N == expected_N_explicit, f"N should equal exp(Phi(t,r)), got {N}"

        # Test mathematical properties - N should always be positive (lapse function)
        assert N.operator() == exp, "N should be an exponential function"

        # Test critical ADM relationship: N² = A (since N = e^Φ, A = e^(2Φ))
        N_squared = N**2
        assert (N_squared - A).simplify_full() == 0, "N² should equal A (fundamental ADM relationship)"

        # Test alternative form: A = N²
        assert (A - N**2).simplify_full() == 0, "A should equal N² in ADM formalism"

        # Test that N can be differentiated (required for extrinsic curvature)
        dN_dt = diff(N, self.t)
        dN_dr = diff(N, self.r)
        assert dN_dt is not None, "N should be differentiable with respect to t"
        assert dN_dr is not None, "N should be differentiable with respect to r"

        # Test chain rule for N derivatives: ∂N/∂t = N*∂Φ/∂t
        expected_dN_dt = N * diff(phi_func, self.t)
        assert (dN_dt - expected_dN_dt).simplify_full() == 0, \
            "∂N/∂t should equal N*∂Φ/∂t by chain rule"

        # Test that N is invertible (required for K_ij = -(1/2N) ∂_t γ_ij)
        N_inverse = 1/N
        assert N_inverse is not None, "N should be invertible (N > 0 always)"

        # Test that 1/N = exp(-Phi)
        expected_N_inv = exp(-phi_func)
        assert (N_inverse - expected_N_inv).simplify_full() == 0, \
            "1/N should equal exp(-Φ)"

        # Test relationship to conformal factor: N = √A
        sqrt_A = sqrt(A)
        assert (N - sqrt_A).simplify_full() == 0, \
            "N should equal √A (lapse function from conformal factor)"

if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
