#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Schwarzschild Equations for GATG
Symbolic representations of Schwarzschild metric equations
"""

from sage.all import *
from typing import Dict, Any, Optional

class schwarzschildEquations:
    """
    Symbolic equations for Schwarzschild metric in GATG
    """

    def __init__(self, vars_dict: Dict[str, Any]):
        """Initialize with shared variable dictionary"""
        self.vars = vars_dict
        self.equations = {}
        self.category_name = 'schwarzschild'

        # Define equations
        self.define_equations()

    def define_equations(self):
        """Define all equations in this category"""
        
        r = self.vars['r']
        M = self.vars['M']
        rs = self.vars['r_s']

        # Symbolic function A(r)
        A = function('A')(r)

        # Schwarzschild ODE
        self.equations['schwarzschild_ode'] = {
            'symbolic': r * diff(A, r) + A - 1,
            'string': "r A'(r) + A(r) - 1 = 0",
            'latex': r"r \frac{dA}{dr} + A - 1 = 0",
            'description': "Schwarzschild vacuum ODE from Einstein equations",
            'solution': 1 - rs/r  # Known solution
        }

        # Schwarzschild constraint (simplified form)
        self.equations['schwarzschild_constraint'] = {
            'symbolic': diff(r * A, r) - 1,
            'string': "(rA)' = 1",
            'latex': r"\frac{d}{dr}(rA) = 1",
            'description': "Schwarzschild constraint equation",
            'solution': 1 - rs/r
        }

        # Schwarzschild solution
        self.equations['schwarzschild_solution'] = {
            'symbolic': 1 - rs/r,
            'string': "A(r) = 1 - r_s/r",
            'latex': r"A(r) = 1 - \frac{r_s}{r}",
            'description': "Schwarzschild metric temporal component solution"
        }

        # Schwarzschild radius relation
        self.equations['schwarzschild_radius'] = {
            'symbolic': 2*self.vars['G']*M/self.vars['c']**2,
            'string': "r_s = 2GM/c²",
            'latex': r"r_s = \frac{2GM}{c^2}",
            'description': "Schwarzschild radius (event horizon for non-rotating black hole)"
        }

        # Schwarzschild metric (full line element)
        self.equations['schwarzschild_metric'] = {
            'symbolic': None,  # Line element is not a single expression
            'string': "ds² = -(1 - r_s/r) dt² + (1 - r_s/r)⁻¹ dr² + r² dΩ²",
            'latex': r"ds^2 = -\left(1 - \frac{r_s}{r}\right) dt^2 + \left(1 - \frac{r_s}{r}\right)^{-1} dr^2 + r^2 d\Omega^2",
            'description': "Schwarzschild metric line element in standard coordinates"
        }

        # Schwarzschild metric (general A(r) form)
        self.equations['schwarzschild_metric_general'] = {
            'symbolic': None,  # Line element is not a single expression
            'string': "ds² = -A(r) dt² + A(r)⁻¹ dr² + r² dΩ²",
            'latex': r"ds^2 = -A(r) dt^2 + A(r)^{-1} dr^2 + r^2 d\Omega^2",
            'description': "Schwarzschild metric in general A(r) form"
        }

        # Alternative Schwarzschild forms
        self.equations['schwarzschild_recovery'] = {
            'symbolic': diff(r * exp(2*function('Phi')(self.vars['t'], r)), r) - 1,
            'string': "(re^(2Φ))' = 1",
            'latex': r"\frac{d}{dr}(re^{2\Phi}) = 1",
            'description': "Schwarzschild recovery equation in exponential form"
        }

        self.equations['schwarzschild_form'] = {
            'symbolic': exp(2*function('Phi')(self.vars['t'], r)) - (1 - 2*M/r),
            'string': "e^(2Φ) = 1 - 2M/r",
            'latex': r"e^{2\Phi} = 1 - \frac{2M}{r}",
            'description': "Schwarzschild solution in exponential lapse form"
        }

        # Vacuum limit condition
        self.equations['vacuum_limit'] = {
            'symbolic': None,  # Limiting condition
            'string': "T^t_r = 0 → ∂_t Φ = 0",
            'latex': r"T^t_r = 0 \to \partial_t \Phi = 0",
            'description': "Vacuum limit: no energy flux implies static temporal potential"
        }



    def verify_schwarzschild_ode(self) -> bool:
        """Verify that the Schwarzschild ODE is satisfied by the known solution"""
        eq_data = self.equations['schwarzschild_ode']
        ode_expr = eq_data['symbolic']
        solution = eq_data['solution']

        # Get variables
        r = self.vars['r']
        A = function('A')(r)

        # Substitute solution into ODE
        # First compute derivative of solution
        A_solution = solution
        dA_dr = diff(A_solution, r)

        # Substitute into ODE: r * A'(r) + A(r) - 1
        result = r * dA_dr + A_solution - 1
        result_simplified = result.simplify_full()

        return result_simplified == 0
    

    def get_equations(self) -> Dict[str, Dict[str, Any]]:
        """Return all equations in this category"""
        return self.equations.copy()

    def get_equation_keys(self) -> list:
        """Return list of equation keys in this category"""
        return list(self.equations.keys())

# Export the class
__all__ = ['schwarzschildEquations']
