#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Cosmological Equations for GATG
Symbolic representations of cosmological (FRW, Friedmann) equations
"""

from sage.all import *
from typing import Dict, Any, Optional

class CosmologicalEquations:
    """
    Symbolic equations for cosmology in GATG
    """

    def __init__(self, vars_dict: Dict[str, Any]):
        """Initialize with shared variable dictionary"""
        self.vars = vars_dict
        self.equations = {}
        self.category_name = 'cosmological'

        # Define equations
        self.define_equations()

    def define_equations(self):
        """Define cosmological equations (FRW, Friedmann, etc.)"""
        t = self.vars['t']
        r = self.vars['r']
        theta = self.vars['theta']
        H = self.vars['H']
        H_0 = self.vars['H_0']
        k = self.vars['k']
        z = self.vars['z']
        rho = self.vars['rho']
        p = self.vars['p']
        G = self.vars['G']
        c = self.vars['c']
        Lambda = self.vars['Lambda']
        rho_m = self.vars['rho_m']
        rho_r = self.vars['rho_r']
        rho_Lambda = self.vars['rho_Lambda']
        Omega_m = self.vars['Omega_m']
        Omega_Lambda = self.vars['Omega_Lambda']

        # Scale factor and temporal potential
        a = function('a')(t)  # Scale factor
        Phi = function('Phi')(t)  # Temporal potential

        # FRW metric equations
        self.equations['frw_metric'] = {
            'symbolic': None,  # Line element is not a single expression
            'string': "ds² = -dt² + a²(t)[dr²/(1-kr²) + r²dΩ²]",
            'latex': r"ds^2 = -dt^2 + a^2(t)\left[\frac{dr^2}{1-kr^2} + r^2 d\Omega^2\right]",
            'description': "Friedmann-Robertson-Walker metric (general curvature k)"
        }

        self.equations['frw_metric_flat'] = {
            'symbolic': None,  # Line element
            'string': "ds² = -dt² + a²(t)[dr² + r²dΩ²]",
            'latex': r"ds^2 = -dt^2 + a^2(t)[dr^2 + r^2 d\Omega^2]",
            'description': "FRW metric for flat universe (k=0)"
        }

        self.equations['frw_lapse_first'] = {
            'symbolic': None,  # Line element
            'string': "ds² = -e^(2Φ)dt² + a²(t)[dr² + r²dΩ²]",
            'latex': r"ds^2 = -e^{2\Phi} dt^2 + a^2(t)[dr^2 + r^2 d\Omega^2]",
            'description': "FRW metric in lapse-first form (GATG approach)"
        }

        # Fundamental relations
        self.equations['hubble_parameter'] = {
            'symbolic': H - diff(a, t)/a,
            'string': "H ≡ ȧ/a",
            'latex': r"H \equiv \frac{\dot{a}}{a}",
            'description': "Hubble parameter definition"
        }

        self.equations['scale_factor_lapse_first'] = {
            'symbolic': a - exp(-Phi),
            'string': "a(t) = e^(-Φ(t))",
            'latex': r"a(t) = e^{-\Phi(t)}",
            'description': "Scale factor in terms of temporal potential (GATG key relation)"
        }

        self.equations['hubble_lapse_first'] = {
            'symbolic': H + diff(Phi, t),
            'string': "H = -∂_t Φ",
            'latex': r"H = -\partial_t \Phi",
            'description': "Hubble parameter in lapse-first formulation"
        }

        # Friedmann equations
        self.equations['friedmann_first'] = {
            'symbolic': H**2 + k/a**2 - (8*pi*G/3)*rho - Lambda/3,
            'string': "H² + k/a² = (8πG/3)ρ + Λ/3",
            'latex': r"H^2 + \frac{k}{a^2} = \frac{8\pi G}{3}\rho + \frac{\Lambda}{3}",
            'description': "First Friedmann equation (energy constraint)"
        }

        self.equations['friedmann_second'] = {
            'symbolic': diff(a, t, 2)/a + (4*pi*G/3)*(rho + 3*p) - Lambda/3,
            'string': "ä/a = -(4πG/3)(ρ + 3p) + Λ/3",
            'latex': r"\frac{\ddot{a}}{a} = -\frac{4\pi G}{3}(\rho + 3p) + \frac{\Lambda}{3}",
            'description': "Second Friedmann equation (acceleration equation)"
        }

        self.equations['friedmann_first_normalized'] = {
            'symbolic': H**2 - (8*pi*G/3)*rho - Lambda/3 + k/a**2,
            'string': "H² = (8πG/3)ρ + Λ/3 - k/a²",
            'latex': r"H^2 = \frac{8\pi G}{3}\rho + \frac{\Lambda}{3} - \frac{k}{a^2}",
            'description': "First Friedmann equation (normalized form)"
        }

        # Matter content equations
        self.equations['matter_dominated'] = {
            'symbolic': rho - rho_m * a**(-3),
            'string': "ρ = ρ_m a^(-3)",
            'latex': r"\rho = \rho_m a^{-3}",
            'description': "Matter density scaling with expansion"
        }

        self.equations['radiation_dominated'] = {
            'symbolic': rho - rho_r * a**(-4),
            'string': "ρ = ρ_r a^(-4)",
            'latex': r"\rho = \rho_r a^{-4}",
            'description': "Radiation density scaling with expansion"
        }

        self.equations['lambda_cdm'] = {
            'symbolic': rho - (rho_m * a**(-3) + rho_Lambda),
            'string': "ρ = ρ_m a^(-3) + ρ_Λ",
            'latex': r"\rho = \rho_m a^{-3} + \rho_\Lambda",
            'description': "ΛCDM matter content (matter + cosmological constant)"
        }

        # Redshift relations
        a_0 = var('a_0', domain='positive')  # Present-day scale factor
        self.equations['redshift_scale_relation'] = {
            'symbolic': (1 + z) - a_0/a,
            'string': "1 + z = a₀/a",
            'latex': r"1 + z = \frac{a_0}{a}",
            'description': "Redshift-scale factor relation"
        }

        # E(z) functions for distance calculations
        E_z = function('E')(z)  # Dimensionless Hubble parameter
        self.equations['ez_lambda_cdm'] = {
            'symbolic': E_z**2 - (Omega_m*(1+z)**3 + Omega_Lambda),
            'string': "E²(z) = Ω_m(1+z)³ + Ω_Λ",
            'latex': r"E^2(z) = \Omega_m(1+z)^3 + \Omega_\Lambda",
            'description': "Dimensionless Hubble parameter for ΛCDM"
        }

        # Distance measures
        d_L = function('d_L')(z)  # Luminosity distance
        d_A = function('d_A')(z)  # Angular diameter distance

        self.equations['luminosity_distance'] = {
            'symbolic': None,  # Involves integral
            'string': "d_L(z) = (1+z) ∫₀^z dz'/E(z')",
            'latex': r"d_L(z) = (1+z) \int_0^z \frac{dz'}{E(z')}",
            'description': "Luminosity distance in terms of redshift"
        }

        self.equations['angular_diameter_distance'] = {
            'symbolic': d_A - d_L/(1+z)**2,
            'string': "d_A(z) = d_L(z)/(1+z)²",
            'latex': r"d_A(z) = \frac{d_L(z)}{(1+z)^2}",
            'description': "Angular diameter distance from luminosity distance"
        }

        # Age of universe
        t_0 = var('t_0', domain='positive')  # Age of universe
        self.equations['age_of_universe'] = {
            'symbolic': None,  # Involves integral
            'string': "t₀ = ∫₀^∞ dz/[(1+z)H₀E(z)]",
            'latex': r"t_0 = \int_0^\infty \frac{dz}{(1+z)H_0 E(z)}",
            'description': "Age of universe from Hubble parameter evolution"
        }

        # Equation of state relations
        w = self.vars['w']
        self.equations['equation_of_state'] = {
            'symbolic': w - p/rho,
            'string': "w = p/ρ",
            'latex': r"w = \frac{p}{\rho}",
            'description': "Equation of state parameter"
        }

        # Critical density
        rho_crit = var('rho_crit', domain='positive')
        self.equations['critical_density'] = {
            'symbolic': rho_crit - 3*H_0**2/(8*pi*G),
            'string': "ρ_crit = 3H₀²/(8πG)",
            'latex': r"\rho_{\text{crit}} = \frac{3H_0^2}{8\pi G}",
            'description': "Critical density for flat universe"
        }

        # Density parameters
        self.equations['omega_matter'] = {
            'symbolic': Omega_m - rho_m/rho_crit,
            'string': "Ω_m = ρ_m/ρ_crit",
            'latex': r"\Omega_m = \frac{\rho_m}{\rho_{\text{crit}}}",
            'description': "Matter density parameter"
        }

        self.equations['omega_lambda'] = {
            'symbolic': Omega_Lambda - rho_Lambda/rho_crit,
            'string': "Ω_Λ = ρ_Λ/ρ_crit",
            'latex': r"\Omega_\Lambda = \frac{\rho_\Lambda}{\rho_{\text{crit}}}",
            'description': "Dark energy density parameter"
        }

        # Flatness condition
        Omega_k = var('Omega_k', domain='real')  # Curvature density parameter
        self.equations['flatness_condition'] = {
            'symbolic': Omega_m + Omega_Lambda + Omega_k - 1,
            'string': "Ω_m + Ω_Λ + Ω_k = 1",
            'latex': r"\Omega_m + \Omega_\Lambda + \Omega_k = 1",
            'description': "Flatness condition (sum of density parameters)"
        }

    def verify_lapse_first_hubble_relation(self) -> bool:
        """Verify consistency between standard and lapse-first Hubble relations"""
        t = self.vars['t']
        a = function('a')(t)
        Phi = function('Phi')(t)
        H = self.vars['H']

        # Standard relation: H = ȧ/a
        H_standard = diff(a, t)/a

        # Lapse-first relation: H = -∂_t Φ
        H_lapse_first = -diff(Phi, t)

        # Scale factor relation: a = e^(-Φ)
        a_from_phi = exp(-Phi)

        # Substitute a = e^(-Φ) into standard relation
        H_from_substitution = diff(a_from_phi, t)/a_from_phi
        H_simplified = H_from_substitution.simplify_full()

        # Should equal -∂_t Φ
        return H_simplified == H_lapse_first

    def verify_matter_scaling(self) -> bool:
        """Verify matter density scaling ρ ∝ a^(-3)"""
        a = function('a')(self.vars['t'])
        rho_m = self.vars['rho_m']

        # Matter density equation: ρ = ρ_m a^(-3)
        matter_density = self.equations['matter_dominated']['symbolic']

        # Check that it reduces to the expected form
        # (This is a structural check)
        return matter_density.has(a**(-3))

    def verify_radiation_scaling(self) -> bool:
        """Verify radiation density scaling ρ ∝ a^(-4)"""
        a = function('a')(self.vars['t'])

        # Radiation density equation: ρ = ρ_r a^(-4)
        radiation_density = self.equations['radiation_dominated']['symbolic']

        # Check that it has the expected a^(-4) dependence
        return radiation_density.has(a**(-4))

    def verify_flatness_condition(self) -> bool:
        """Verify that density parameters sum to unity for flat universe"""
        Omega_m = self.vars['Omega_m']
        Omega_Lambda = self.vars['Omega_Lambda']

        # Get flatness condition: Ω_m + Ω_Λ + Ω_k = 1
        flatness = self.equations['flatness_condition']['symbolic']

        # For flat universe (k=0), we should have Ω_k = 0
        # So Ω_m + Ω_Λ = 1
        flat_condition = flatness.substitute({var('Omega_k'): 0})
        expected = Omega_m + Omega_Lambda - 1

        return flat_condition == expected

    def verify_cosmological_redshift_relation(self) -> bool:
        """Verify redshift-scale factor relation"""
        z = self.vars['z']
        a = function('a')(self.vars['t'])
        a_0 = var('a_0', domain='positive')

        # Get redshift relation: 1 + z = a₀/a
        redshift_rel = self.equations['redshift_scale_relation']['symbolic']

        # Check that it has the expected structure
        expected = (1 + z) - a_0/a

        return redshift_rel == expected

    def get_equations(self) -> Dict[str, Dict[str, Any]]:
        """Return all equations in this category"""
        return self.equations.copy()

    def get_equation_keys(self) -> list:
        """Return list of equation keys in this category"""
        return list(self.equations.keys())

# Export the class
__all__ = ['CosmologicalEquations']