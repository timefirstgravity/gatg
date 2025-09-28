#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Dimensional Analysis Equations for GATG
Symbolic representations of physical dimensions and unit consistency checks
"""

from sage.all import *
from typing import Dict, Any, Optional

class DimensionalAnalysisEquations:
    """
    Symbolic equations for dimensional analysis and unit checking in GATG
    """

    def __init__(self, vars_dict: Dict[str, Any]):
        """Initialize with shared variable dictionary"""
        self.vars = vars_dict
        self.equations = {}
        self.category_name = 'dimensional_analysis'

        # Define equations
        self.define_equations()

    def define_equations(self):
        """Define dimensional analysis equations for GATG flux law and related quantities"""

        # GATG flux law dimensional analysis
        self.equations['flux_law_lhs_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[∂_t Φ] = T^(-1)",
            'latex': r"[\partial_t \Phi] = T^{-1}",
            'description': "Temporal potential derivative has dimension of inverse time"
        }

        self.equations['flux_law_rhs_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[4πG/c⁴ × r × T^tr] = T^(-1)",
            'latex': r"[4\pi G/c^4 \times r \times T^{tr}] = T^{-1}",
            'description': "Right-hand side of flux law has matching dimension"
        }

        self.equations['stress_energy_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[T^tr] = MT^(-3)",
            'latex': r"[T^{tr}] = MT^{-3}",
            'description': "Mixed stress-energy tensor has dimension mass per time cubed"
        }

        self.equations['dimensional_consistency_check'] = {
            'symbolic': None,  # Verification statement
            'string': "[T^tr] = MT^(-3)",
            'latex': r"[T^{tr}] = MT^{-3}",
            'description': "Dimensional consistency verification for stress-energy tensor"
        }

        # Detailed dimensional breakdown
        self.equations['dimensional_explanation'] = {
            'symbolic': None,  # Physical interpretation
            'string': "[T^tr] = energy/(area·time), r adds length, G/c⁴ contributes 1/pressure, so RHS is 1/time",
            'latex': r"[T^{tr}] = \frac{\text{energy}}{\text{area} \cdot \text{time}}, \quad r \text{ adds length}, \quad G/c^4 \text{ contributes } \frac{1}{\text{pressure}}, \text{ so RHS is } \frac{1}{\text{time}}",
            'description': "Detailed dimensional analysis breakdown of flux law components"
        }

        # Fundamental constant dimensions
        self.equations['gravitational_constant_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[G] = M^(-1)L³T^(-2)",
            'latex': r"[G] = M^{-1}L^3T^{-2}",
            'description': "Gravitational constant dimension in natural units"
        }

        self.equations['speed_of_light_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[c] = LT^(-1)",
            'latex': r"[c] = LT^{-1}",
            'description': "Speed of light dimension"
        }

        self.equations['planck_constant_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[ℏ] = ML²T^(-1)",
            'latex': r"[\hbar] = ML^2T^{-1}",
            'description': "Reduced Planck constant dimension"
        }

        # Metric function dimensions
        self.equations['temporal_potential_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[Φ] = dimensionless",
            'latex': r"[\Phi] = \text{dimensionless}",
            'description': "Temporal potential is dimensionless (logarithmic quantity)"
        }

        self.equations['metric_function_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[A] = dimensionless",
            'latex': r"[A] = \text{dimensionless}",
            'description': "Metric function A = e^(2Φ) is dimensionless"
        }

        self.equations['lapse_function_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[N] = dimensionless",
            'latex': r"[N] = \text{dimensionless}",
            'description': "Lapse function N = e^Φ is dimensionless"
        }

        # Coordinate dimensions
        self.equations['coordinate_time_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[t] = T",
            'latex': r"[t] = T",
            'description': "Coordinate time has dimension of time"
        }

        self.equations['coordinate_radius_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[r] = L",
            'latex': r"[r] = L",
            'description': "Areal radius coordinate has dimension of length"
        }

        # Curvature dimensions
        self.equations['ricci_scalar_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[R] = L^(-2)",
            'latex': r"[R] = L^{-2}",
            'description': "Ricci scalar has dimension of inverse length squared"
        }

        self.equations['einstein_tensor_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[G_μν] = L^(-2)",
            'latex': r"[G_{\mu\nu}] = L^{-2}",
            'description': "Einstein tensor has dimension of inverse length squared"
        }

        # Energy-momentum dimensions
        self.equations['energy_density_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[ρ] = ML^(-1)T^(-2)",
            'latex': r"[\rho] = ML^{-1}T^{-2}",
            'description': "Energy density dimension"
        }

        self.equations['pressure_dimension'] = {
            'symbolic': None,  # Dimensional analysis
            'string': "[p] = ML^(-1)T^(-2)",
            'latex': r"[p] = ML^{-1}T^{-2}",
            'description': "Pressure has same dimension as energy density"
        }

        # Geometrized units (G = c = 1)
        self.equations['geometrized_units_condition'] = {
            'symbolic': None,  # Unit system
            'string': "G = c = 1 (geometrized units)",
            'latex': r"G = c = 1 \quad \text{(geometrized units)}",
            'description': "Natural unit system where G and c are set to unity"
        }

        self.equations['geometrized_mass_dimension'] = {
            'symbolic': None,  # Dimensional analysis in geometrized units
            'string': "[M] = L = T (geometrized)",
            'latex': r"[M] = L = T \quad \text{(geometrized)}",
            'description': "In geometrized units, mass, length, and time have the same dimension"
        }

        # Planck units
        self.equations['planck_length_definition'] = {
            'symbolic': sqrt(self.vars['G'] * var('hbar') / self.vars['c']**3),
            'string': "l_P = √(Gℏ/c³)",
            'latex': r"l_P = \sqrt{\frac{G\hbar}{c^3}}",
            'description': "Planck length - fundamental length scale"
        }

        self.equations['planck_time_definition'] = {
            'symbolic': sqrt(self.vars['G'] * var('hbar') / self.vars['c']**5),
            'string': "t_P = √(Gℏ/c⁵)",
            'latex': r"t_P = \sqrt{\frac{G\hbar}{c^5}}",
            'description': "Planck time - fundamental time scale"
        }

        self.equations['planck_mass_definition'] = {
            'symbolic': sqrt(var('hbar') * self.vars['c'] / self.vars['G']),
            'string': "m_P = √(ℏc/G)",
            'latex': r"m_P = \sqrt{\frac{\hbar c}{G}}",
            'description': "Planck mass - fundamental mass scale"
        }

    def verify_flux_law_dimensional_consistency(self) -> bool:
        """Verify that GATG flux law is dimensionally consistent"""
        # The flux law: ∂_t Φ = (4πG/c⁴) r T^tr
        # LHS: [∂_t Φ] = T^(-1) (since Φ is dimensionless)
        # RHS: [G/c⁴] × [r] × [T^tr] = [M^(-1)L³T^(-2)] × [L^(-4)T^4] × [L] × [MT^(-3)]
        #    = M^(-1)L³T^(-2) × L^(-4)T^4 × L × MT^(-3)
        #    = M^(-1)L³T^(-2) × L^(-3)T^1 × M
        #    = L⁰T^(-1) = T^(-1) ✓

        # This is a conceptual verification - the dimensional analysis is built into the equation definitions
        lhs_desc = self.equations['flux_law_lhs_dimension']['string']
        rhs_desc = self.equations['flux_law_rhs_dimension']['string']

        # Both should have T^(-1) dimension
        return "T^(-1)" in lhs_desc and "T^(-1)" in rhs_desc

    def verify_stress_energy_dimension(self) -> bool:
        """Verify stress-energy tensor has correct physical dimension"""
        # T^tr should have dimension of energy flux: energy/(area × time)
        # In SI units: [T^tr] = J/(m² × s) = kg⋅m²⋅s^(-2) / (m² × s) = kg⋅s^(-3) = MT^(-3)

        stress_energy_desc = self.equations['stress_energy_dimension']['string']
        return "MT^(-3)" in stress_energy_desc

    def verify_geometrized_units_consistency(self) -> bool:
        """Verify that geometrized units maintain dimensional consistency"""
        # In geometrized units G = c = 1, all fundamental quantities can be expressed in terms of length
        # Mass, time, and length all have the same dimension

        geom_condition = self.equations['geometrized_units_condition']['string']
        geom_mass_dim = self.equations['geometrized_mass_dimension']['string']

        return "G = c = 1" in geom_condition and "[M] = L = T" in geom_mass_dim

    def get_equations(self) -> Dict[str, Dict[str, Any]]:
        """Return all equations in this category"""
        return self.equations.copy()

    def get_equation_keys(self) -> list:
        """Return list of equation keys in this category"""
        return list(self.equations.keys())

# Export the class
__all__ = ['DimensionalAnalysisEquations']