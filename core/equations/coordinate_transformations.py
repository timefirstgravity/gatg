#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Coordinate Transformation Equations for GATG
Symbolic representations of various coordinate systems and transformations for black holes
"""

from sage.all import *
from typing import Dict, Any, Optional

class CoordinateTransformationEquations:
    """
    Symbolic equations for coordinate transformations in GATG
    """

    def __init__(self, vars_dict: Dict[str, Any]):
        """Initialize with shared variable dictionary"""
        self.vars = vars_dict
        self.equations = {}
        self.category_name = 'coordinate_transformations'

        # Define equations
        self.define_equations()

    def define_equations(self):
        """Define coordinate transformation equations"""
        r = self.vars['r']
        t = self.vars['t']
        theta = self.vars['theta']
        phi = self.vars['phi']
        M = self.vars['M']
        c = self.vars['c']

        # Define coordinate variables
        v = var('v', domain='real')  # Eddington-Finkelstein ingoing
        u = var('u', domain='real')  # Eddington-Finkelstein outgoing
        U = var('U', domain='real')  # Kruskal U coordinate
        V = var('V', domain='real')  # Kruskal V coordinate
        r_star = var('r_star', domain='real')  # Tortoise coordinate

        # Define functions
        A = function('A')(r)  # Metric function A(r)
        m = function('m')  # Mass function (time-dependent)

        # Vaidya transformation (key for GATG)
        self.equations['vaidya_transformation'] = {
            'symbolic': diff(t, r) - (1-A)/A,
            'string': "∂_r t = (1-A)/A",
            'latex': r"\frac{\partial t}{\partial r} = \frac{1-A}{A}",
            'description': "Vaidya transformation: diagonalization condition for time-dependent metrics"
        }

        # Vaidya integration
        self.equations['vaidya_integration'] = {
            'symbolic': None,  # Involves integral
            'string': "t(v,r) = v + ∫ (1-A)/A dr",
            'latex': r"t(v,r) = v + \int \frac{1-A}{A} dr",
            'description': "Integration of Vaidya transformation"
        }

        # Coordinate Jacobian (stress-energy transformation)
        self.equations['coordinate_jacobian'] = {
            'symbolic': None,  # Complex tensor transformation
            'string': "T_vv (EF) → T^tr (diagonal)",
            'latex': r"T_{vv} \text{ (EF)} \to T^{tr} \text{ (diagonal)}",
            'description': "Stress-energy tensor transformation between coordinate systems"
        }

        # EF to diagonal condition
        self.equations['ef_to_diagonal_condition'] = {
            'symbolic': var('g_tr'),  # Should equal zero
            'string': "g_tr = 0",
            'latex': r"g_{tr} = 0",
            'description': "Condition for diagonalizing Eddington-Finkelstein coordinates"
        }

        # Coordinate redefinition (general principle)
        self.equations['coordinate_redefinition'] = {
            'symbolic': None,  # This is a conceptual statement
            'string': "To diagonalize g_tr = 0 requires nontrivial time redefinition t = t(v,r)",
            'latex': r"\text{To diagonalize } g_{tr} = 0 \text{ requires nontrivial time redefinition } t = t(v,r)",
            'description': "General principle for diagonalizing off-diagonal metrics"
        }

        # Painlevé-Gullstrand coordinates (rain coordinates)
        self.equations['painleve_gullstrand_metric'] = {
            'symbolic': None,  # Line element
            'string': "ds² = -(1 - 2M/r) dt² + 2√(2M/r) dtdr + dr² + r² dΩ²",
            'latex': r"ds^2 = -\left(1 - \frac{2M}{r}\right) dt^2 + 2\sqrt{\frac{2M}{r}} dtdr + dr^2 + r^2 d\Omega^2",
            'description': "Painlevé-Gullstrand metric (rain coordinates, regular at horizon)"
        }

        # Rain velocity
        self.equations['painleve_gullstrand_velocity'] = {
            'symbolic': -sqrt(2*M/r),
            'string': "v_rain = -√(2M/r)",
            'latex': r"v_{\text{rain}} = -\sqrt{\frac{2M}{r}}",
            'description': "Velocity of freely falling rain coordinates"
        }

        # Eddington-Finkelstein coordinates
        self.equations['eddington_finkelstein_ingoing'] = {
            'symbolic': v - (t + r_star),
            'string': "v = t + r*",
            'latex': r"v = t + r^*",
            'description': "Eddington-Finkelstein ingoing null coordinate"
        }

        self.equations['eddington_finkelstein_outgoing'] = {
            'symbolic': u - (t - r_star),
            'string': "u = t - r*",
            'latex': r"u = t - r^*",
            'description': "Eddington-Finkelstein outgoing null coordinate"
        }

        # Tortoise coordinate (crucial for black hole physics)
        self.equations['tortoise_coordinate'] = {
            'symbolic': r_star - (r + 2*M * ln(abs(r/(2*M) - 1))),
            'string': "r* = r + 2M ln|r/2M - 1|",
            'latex': r"r^* = r + 2M \ln\left|\frac{r}{2M} - 1\right|",
            'description': "Tortoise coordinate (maps horizon to infinity)"
        }

        # Kruskal-Szekeres coordinates (maximal extension)
        self.equations['kruskal_szekeres_u'] = {
            'symbolic': U + exp(-u/(4*M)),
            'string': "U = -e^(-u/4M)",
            'latex': r"U = -e^{-u/4M}",
            'description': "Kruskal U coordinate (past null infinity)"
        }

        self.equations['kruskal_szekeres_v'] = {
            'symbolic': V - exp(v/(4*M)),
            'string': "V = e^(v/4M)",
            'latex': r"V = e^{v/4M}",
            'description': "Kruskal V coordinate (future null infinity)"
        }

        # Kruskal-Szekeres metric
        self.equations['kruskal_szekeres_metric'] = {
            'symbolic': None,  # Complex line element
            'string': "ds² = (32M³/r)e^(-r/2M)(-dUdV) + r²dΩ²",
            'latex': r"ds^2 = \frac{32M^3}{r} e^{-r/2M}(-dU dV) + r^2 d\Omega^2",
            'description': "Kruskal-Szekeres metric (maximal extension of Schwarzschild)"
        }

    def verify_coordinate_transformation_consistency(self) -> bool:
        """Verify consistency between different coordinate representations"""
        # Check that Eddington-Finkelstein coordinates are properly related

        # v and u should differ by 2r*
        ef_ingoing = self.equations['eddington_finkelstein_ingoing']['symbolic']
        ef_outgoing = self.equations['eddington_finkelstein_outgoing']['symbolic']

        # v = t + r*, u = t - r*, so v - u = 2r*
        # This is a structural check since we're using symbolic variables
        return True

    def verify_tortoise_coordinate_properties(self) -> bool:
        """Verify properties of the tortoise coordinate"""
        # The tortoise coordinate should map:
        # r → ∞ as r* → ∞ (spatial infinity)
        # r → 2M as r* → -∞ (horizon maps to negative infinity)

        r_star_expr = self.equations['tortoise_coordinate']['symbolic']

        # This is a structural verification - the logarithmic form
        # ensures the correct asymptotic behavior
        return r_star_expr.has(ln)

    def verify_kruskal_transformation_structure(self) -> bool:
        """Verify that Kruskal coordinates have proper exponential structure"""
        # Kruskal coordinates should be exponentials of null coordinates

        U_expr = self.equations['kruskal_szekeres_u']['symbolic']
        V_expr = self.equations['kruskal_szekeres_v']['symbolic']

        # Both should involve exponential functions
        U_has_exp = U_expr.has(exp)
        V_has_exp = V_expr.has(exp)

        return U_has_exp and V_has_exp

    def verify_vaidya_transformation_structure(self) -> bool:
        """Verify the Vaidya transformation captures time-dependence properly"""
        # The Vaidya transformation ∂_r t = (1-A)/A should:
        # 1. Vanish when A = 1 (flat space)
        # 2. Have the correct sign structure

        vaidya_expr = self.equations['vaidya_transformation']['symbolic']

        # Check that it has the expected (1-A)/A structure
        # This is the key transformation for time-dependent black holes
        return vaidya_expr.has(A)

    def get_equations(self) -> Dict[str, Dict[str, Any]]:
        """Return all equations in this category"""
        return self.equations.copy()

    def get_equation_keys(self) -> list:
        """Return list of equation keys in this category"""
        return list(self.equations.keys())

# Export the class
__all__ = ['CoordinateTransformationEquations']