#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Flux Law Equations for GATG
Symbolic representations of flux law equations
"""

from sage.all import *
from typing import Dict, Any, Optional

class FluxLawEquations:
    """
    Symbolic equations for flux law in GATG
    """

    def __init__(self, vars_dict: Dict[str, Any]):
        """Initialize with shared variable dictionary"""
        self.vars = vars_dict
        self.equations = {}
        self.category_name = 'flux_law'

        # Define equations
        self.define_equations()

    def define_equations(self):
        """Define all equations in this category"""
        
        r = self.vars['r']
        t = self.vars['t']
        G = self.vars['G']
        c = self.vars['c']

        # Define symbolic functions
        Phi = function('Phi')(t, r)
        T_tr = var('T_tr')  # Stress-energy component

        # Primary flux law
        self.equations['flux_law'] = {
            'symbolic': diff(Phi, t) == (4*pi*G/c**4) * r * T_tr,
            'string': "∂_t Φ = (4πG/c⁴) r T^t_r",
            'latex': r"\partial_t \Phi = \frac{4\pi G}{c^4} r T^t_r",
            'description': "Temporal flux law: time evolution driven by energy flux"
        }

        # Flux law in geometrized units (G=c=1)
        self.equations['flux_law_geometrized'] = {
            'symbolic': diff(Phi, t) == 4*pi * r * T_tr,
            'string': "∂_t Φ = 4π r T^t_r",
            'latex': r"\partial_t \Phi = 4\pi r T^t_r",
            'description': "Flux law in geometrized units (G=c=1)"
        }

        # Einstein tensor component relation
        G_tr = var('G_tr')
        self.equations['einstein_tr'] = {
            'symbolic': G_tr == (2/r) * diff(Phi, t),
            'string': "G^t_r = (2/r) ∂_t Φ",
            'latex': r"G^t_r = \frac{2}{r} \partial_t \Phi",
            'description': "Einstein tensor mixed component in terms of temporal flux"
        }

        # Einstein field equation constraint
        self.equations['einstein_constraint'] = {
            'symbolic': G_tr == (8*pi*G/c**4) * T_tr,
            'string': "G^t_r = (8πG/c⁴) T^t_r",
            'latex': r"G^t_r = \frac{8\pi G}{c^4} T^t_r",
            'description': "Einstein field equation for mixed time-radial component"
        }

        # Flux law with mixed indices (includes metric factor)
        self.equations['flux_law_mixed_indices'] = {
            'symbolic': diff(Phi, t) == (4*pi*G/c**4) * r * exp(2*Phi) * T_tr,
            'string': "∂_t Φ = (4πG/c⁴) r e^(2Φ) T^t_r",
            'latex': r"\partial_t \Phi = \frac{4\pi G}{c^4} r e^{2\Phi} T^t_r",
            'description': "Flux law using mixed index stress-energy tensor with metric factor"
        }
        


    def verify_flux_law_consistency(self) -> bool:
        """Verify consistency between flux law and Einstein tensor relations"""
        # Verify: If ∂_t Φ = (4πG/c⁴) r T^t_r and G^t_r = (2/r) ∂_t Φ
        # Then G^t_r = (2/r) × (4πG/c⁴) r T^t_r = (8πG/c⁴) T^t_r
        # This should match the Einstein constraint G^t_r = (8πG/c⁴) T^t_r

        r = self.vars['r']
        G = self.vars['G']
        c = self.vars['c']
        T_tr = var('T_tr')

        # Direct calculation from flux law
        flux_term = (4*pi*G/c**4) * r * T_tr
        G_tr_from_flux = (2/r) * flux_term

        # Expected from Einstein constraint
        G_tr_expected = (8*pi*G/c**4) * T_tr

        # Check consistency
        difference = (G_tr_from_flux - G_tr_expected).simplify_full()

        return difference == 0
    

    def get_equations(self) -> Dict[str, Dict[str, Any]]:
        """Return all equations in this category"""
        return self.equations.copy()

    def get_equation_keys(self) -> list:
        """Return list of equation keys in this category"""
        return list(self.equations.keys())

# Export the class
__all__ = ['FluxLawEquations']
