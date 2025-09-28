#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Kerr-Newman Solution Equations for GATG
Symbolic representations of charged rotating black hole equations
"""

from sage.all import *
from typing import Dict, Any, Optional

class KerrNewmanEquations:
    """
    Symbolic equations for Kerr-Newman solutions (charged rotating black holes) in GATG
    """

    def __init__(self, vars_dict: Dict[str, Any]):
        """Initialize with shared variable dictionary"""
        self.vars = vars_dict
        self.equations = {}
        self.category_name = 'kerr_newman'

        # Define equations
        self.define_equations()

    def define_equations(self):
        """Define Kerr-Newman solution equations (charged rotating black hole)"""
        r = self.vars['r']
        theta = self.vars['theta']
        M = self.vars['M']
        a = self.vars['a']  # Kerr parameter (angular momentum)
        Q = self.vars['Q']  # Electric charge

        # Kerr-Newman Δ function
        self.equations['kerr_newman_delta'] = {
            'symbolic': r**2 - 2*M*r + a**2 + Q**2,
            'string': "Δ = r² - 2Mr + a² + Q²",
            'latex': r"\Delta = r^2 - 2Mr + a^2 + Q^2",
            'description': "Kerr-Newman Δ function (includes both rotation and charge)"
        }

        # Event horizons
        self.equations['kerr_newman_horizons'] = {
            'symbolic': [M + sqrt(M**2 - a**2 - Q**2), M - sqrt(M**2 - a**2 - Q**2)],
            'string': "r_± = M ± √(M² - a² - Q²)",
            'latex': r"r_\pm = M \pm \sqrt{M^2 - a^2 - Q^2}",
            'description': "Kerr-Newman event horizons (rotation and charge effects)"
        }

        # Σ function (same as Kerr)
        self.equations['kerr_newman_sigma'] = {
            'symbolic': r**2 + a**2 * cos(theta)**2,
            'string': "Σ = r² + a²cos²θ",
            'latex': r"\Sigma = r^2 + a^2\cos^2\theta",
            'description': "Kerr-Newman Σ function (same as pure Kerr)"
        }

        # Existence condition (no naked singularities)
        self.equations['kerr_newman_existence_condition'] = {
            'symbolic': M**2 - a**2 - Q**2,
            'string': "M² ≥ a² + Q² (no naked singularity)",
            'latex': r"M^2 \geq a^2 + Q^2 \quad \text{(no naked singularity)}",
            'description': "Physical constraint preventing naked singularities in Kerr-Newman"
        }

        # Extremal condition
        self.equations['kerr_newman_extremal'] = {
            'symbolic': M**2 - a**2 - Q**2,
            'string': "M² = a² + Q² (extremal case)",
            'latex': r"M^2 = a^2 + Q^2 \quad \text{(extremal case)}",
            'description': "Extremal Kerr-Newman: horizons coincide"
        }

        # Limits and reductions
        self.equations['kerr_newman_to_kerr_limit'] = {
            'symbolic': None,  # Conceptual limit
            'string': "Q → 0: Kerr-Newman → Kerr",
            'latex': r"Q \to 0: \text{Kerr-Newman} \to \text{Kerr}",
            'description': "Uncharged limit reduces to pure Kerr solution"
        }

        self.equations['kerr_newman_to_reissner_nordstrom_limit'] = {
            'symbolic': None,  # Conceptual limit
            'string': "a → 0: Kerr-Newman → Reissner-Nordström",
            'latex': r"a \to 0: \text{Kerr-Newman} \to \text{Reissner-Nordström}",
            'description': "Non-rotating limit reduces to charged Reissner-Nordström"
        }

        self.equations['kerr_newman_to_schwarzschild_limit'] = {
            'symbolic': None,  # Conceptual limit
            'string': "a → 0, Q → 0: Kerr-Newman → Schwarzschild",
            'latex': r"a \to 0, Q \to 0: \text{Kerr-Newman} \to \text{Schwarzschild}",
            'description': "Non-rotating, uncharged limit reduces to Schwarzschild"
        }

        # Physical interpretation
        self.equations['kerr_newman_parameters'] = {
            'symbolic': None,  # Parameter description
            'string': "M = mass, a = J/(Mc) = rotation, Q = charge",
            'latex': r"M = \text{mass}, \quad a = \frac{J}{Mc} = \text{rotation}, \quad Q = \text{charge}",
            'description': "Three fundamental parameters of Kerr-Newman black holes"
        }

        # Ergosphere (same formula as Kerr since charge doesn't affect ergosphere)
        self.equations['kerr_newman_ergosphere'] = {
            'symbolic': M + sqrt(M**2 - a**2 * cos(theta)**2),
            'string': "r_e = M + √(M² - a²cos²θ)",
            'latex': r"r_e = M + \sqrt{M^2 - a^2\cos^2\theta}",
            'description': "Kerr-Newman ergosphere (charge doesn't affect ergosphere boundary)"
        }

    def verify_kerr_newman_horizons(self) -> bool:
        """Verify that horizons satisfy Δ = 0"""
        M = self.vars['M']
        a = self.vars['a']
        Q = self.vars['Q']
        r = self.vars['r']

        # Get Δ function and horizon locations
        delta_function = self.equations['kerr_newman_delta']['symbolic']
        horizons = self.equations['kerr_newman_horizons']['symbolic']

        r_plus = horizons[0]   # M + sqrt(M^2 - a^2 - Q^2)
        r_minus = horizons[1]  # M - sqrt(M^2 - a^2 - Q^2)

        # Check that Δ(r_±) = 0
        delta_at_r_plus = delta_function.substitute({r: r_plus}).simplify_full()
        delta_at_r_minus = delta_function.substitute({r: r_minus}).simplify_full()

        return delta_at_r_plus == 0 and delta_at_r_minus == 0

    def verify_kerr_newman_limits(self) -> bool:
        """Verify that Kerr-Newman reduces to known solutions in appropriate limits"""
        M = self.vars['M']
        a = self.vars['a']
        Q = self.vars['Q']

        # Get Kerr-Newman Δ function
        delta_kn = self.equations['kerr_newman_delta']['symbolic']

        # Q → 0 limit should give Kerr: Δ = r² - 2Mr + a²
        delta_kerr = delta_kn.substitute({Q: 0}).simplify_full()
        expected_kerr = self.vars['r']**2 - 2*M*self.vars['r'] + a**2

        # a → 0 limit should give Reissner-Nordström: Δ = r² - 2Mr + Q²
        delta_rn = delta_kn.substitute({a: 0}).simplify_full()
        expected_rn = self.vars['r']**2 - 2*M*self.vars['r'] + Q**2

        # Both limits → Schwarzschild: Δ = r² - 2Mr
        delta_schwarzschild = delta_kn.substitute({a: 0, Q: 0}).simplify_full()
        expected_schwarzschild = self.vars['r']**2 - 2*M*self.vars['r']

        return (delta_kerr == expected_kerr and
                delta_rn == expected_rn and
                delta_schwarzschild == expected_schwarzschild)

    def verify_extremal_condition(self) -> bool:
        """Verify properties of extremal Kerr-Newman black holes"""
        M = self.vars['M']
        a = self.vars['a']
        Q = self.vars['Q']

        # In extremal case, M² = a² + Q²
        # Both horizons should coincide at r = M

        extremal_condition = self.equations['kerr_newman_extremal']['symbolic']
        horizons = self.equations['kerr_newman_horizons']['symbolic']

        # Under extremal condition M² - a² - Q² = 0, both horizons are at r = M
        r_plus_extremal = horizons[0].substitute({M**2 - a**2 - Q**2: 0}).simplify_full()
        r_minus_extremal = horizons[1].substitute({M**2 - a**2 - Q**2: 0}).simplify_full()

        return r_plus_extremal == M and r_minus_extremal == M

    def get_equations(self) -> Dict[str, Dict[str, Any]]:
        """Return all equations in this category"""
        return self.equations.copy()

    def get_equation_keys(self) -> list:
        """Return list of equation keys in this category"""
        return list(self.equations.keys())

# Export the class
__all__ = ['KerrNewmanEquations']