#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Reissner-Nordström Equations for GATG
Symbolic representations of charged black hole equations
"""

from sage.all import *
from typing import Dict, Any, Optional

class ReissnerNordstromEquations:
    """
    Symbolic equations for Reissner-Nordström (charged) black holes in GATG
    """

    def __init__(self, vars_dict: Dict[str, Any]):
        """Initialize with shared variable dictionary"""
        self.vars = vars_dict
        self.equations = {}
        self.category_name = 'reissner_nordstrom'

        # Define equations
        self.define_equations()

    def define_equations(self):
        """Define Reissner-Nordström metric equations (charged black hole)"""
        r = self.vars['r']
        theta = self.vars['theta']
        M = self.vars['M']
        Q = self.vars['Q']  # Electric charge

        # Reissner-Nordström metric
        self.equations['reissner_nordstrom_metric'] = {
            'symbolic': None,  # Line element
            'string': "ds² = -(1 - 2M/r + Q²/r²) dt² + (1 - 2M/r + Q²/r²)⁻¹ dr² + r² dΩ²",
            'latex': r"ds^2 = -\left(1 - \frac{2M}{r} + \frac{Q^2}{r^2}\right) dt^2 + \left(1 - \frac{2M}{r} + \frac{Q^2}{r^2}\right)^{-1} dr^2 + r^2 d\Omega^2",
            'description': "Reissner-Nordström metric for charged black hole"
        }

        # Metric function A(r) for charged black hole
        self.equations['reissner_nordstrom_solution'] = {
            'symbolic': 1 - 2*M/r + Q**2/r**2,
            'string': "A(r) = 1 - 2M/r + Q²/r²",
            'latex': r"A(r) = 1 - \frac{2M}{r} + \frac{Q^2}{r^2}",
            'description': "Temporal metric function for Reissner-Nordström black hole"
        }

        # Event horizons (where A(r) = 0)
        self.equations['reissner_nordstrom_horizons'] = {
            'symbolic': [M + sqrt(M**2 - Q**2), M - sqrt(M**2 - Q**2)],
            'string': "r_± = M ± √(M² - Q²)",
            'latex': r"r_\pm = M \pm \sqrt{M^2 - Q^2}",
            'description': "Inner and outer horizons of Reissner-Nordström black hole"
        }

        # Electric charge (parameter)
        self.equations['reissner_nordstrom_charge'] = {
            'symbolic': Q,
            'string': "Q = electric charge",
            'latex': r"Q = \text{electric charge}",
            'description': "Electric charge parameter of the black hole"
        }

        # Extremal condition
        self.equations['reissner_nordstrom_extremal'] = {
            'symbolic': M**2 - Q**2,
            'string': "M² - Q² = 0 (extremal case)",
            'latex': r"M^2 - Q^2 = 0 \quad \text{(extremal case)}",
            'description': "Extremal condition where inner and outer horizons coincide"
        }

        # Physical constraint (no naked singularities)
        self.equations['reissner_nordstrom_constraint'] = {
            'symbolic': M**2 - Q**2,
            'string': "M² ≥ Q² (cosmic censorship)",
            'latex': r"M^2 \geq Q^2 \quad \text{(cosmic censorship)}",
            'description': "Physical constraint preventing naked singularities"
        }

    def verify_reissner_nordstrom_horizons(self) -> bool:
        """Verify that horizons satisfy A(r) = 0"""
        M = self.vars['M']
        Q = self.vars['Q']
        r = self.vars['r']

        # Get metric function and horizon locations
        A_function = self.equations['reissner_nordstrom_solution']['symbolic']
        horizons = self.equations['reissner_nordstrom_horizons']['symbolic']

        r_plus = horizons[0]   # M + sqrt(M^2 - Q^2)
        r_minus = horizons[1]  # M - sqrt(M^2 - Q^2)

        # Check that A(r_±) = 0
        A_at_r_plus = A_function.substitute({r: r_plus}).simplify_full()
        A_at_r_minus = A_function.substitute({r: r_minus}).simplify_full()

        return A_at_r_plus == 0 and A_at_r_minus == 0

    def verify_schwarzschild_limit(self) -> bool:
        """Verify that Q → 0 gives Schwarzschild solution"""
        M = self.vars['M']
        Q = self.vars['Q']
        r = self.vars['r']

        # Get Reissner-Nordström solution
        A_rn = self.equations['reissner_nordstrom_solution']['symbolic']

        # Take Q → 0 limit
        A_schwarzschild = A_rn.substitute({Q: 0}).simplify_full()

        # Should equal 1 - 2M/r
        expected = 1 - 2*M/r

        return A_schwarzschild == expected

    def verify_extremal_condition(self) -> bool:
        """Verify properties of extremal black holes"""
        M = self.vars['M']
        Q = self.vars['Q']

        # In extremal case, M² = Q²
        # Both horizons should coincide at r = M

        extremal_condition = self.equations['reissner_nordstrom_extremal']['symbolic']

        # When M² - Q² = 0, both horizons are at r = M
        horizons = self.equations['reissner_nordstrom_horizons']['symbolic']

        # Under extremal condition, r_+ = r_- = M
        r_plus_extremal = horizons[0].substitute({Q**2: M**2}).simplify_full()
        r_minus_extremal = horizons[1].substitute({Q**2: M**2}).simplify_full()

        return r_plus_extremal == M and r_minus_extremal == M

    def get_equations(self) -> Dict[str, Dict[str, Any]]:
        """Return all equations in this category"""
        return self.equations.copy()

    def get_equation_keys(self) -> list:
        """Return list of equation keys in this category"""
        return list(self.equations.keys())

# Export the class
__all__ = ['ReissnerNordstromEquations']