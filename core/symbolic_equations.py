#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Symbolic Equations - Modular Architecture
Connects string equation documentation with actual symbolic mathematics

This module provides a rigorous connection between the documented equations
in GATGEquations and their symbolic SageMath representations, enabling:
- Mathematical verification of documented equations
- Computation using symbolic expressions
- Automatic generation of LaTeX from symbolic forms
- Verification of mathematical identities and relationships

The equations are now organized into category-specific modules for better
maintainability and organization.
"""

from sage.all import *
from typing import Dict, Any, Optional, Tuple
import os
import sys

# Import all category modules
try:
    from .equations.schwarzschild import schwarzschildEquations
    from .equations.flux_law import FluxLawEquations
    from .equations.adm import AdmEquations
    from .equations.kerr_equations import KerrEquations
    from .equations.cosmological import CosmologicalEquations
    from .equations.christoffel import ChristoffelEquations
    from .equations.conservation import ConservationEquations
    from .equations.coordinate_transformations import CoordinateTransformationEquations
    from .equations.reissner_nordstrom import ReissnerNordstromEquations
    from .equations.vaidya import VaidyaEquations
    from .equations.kerr_newman import KerrNewmanEquations
    from .equations.dimensional_analysis import DimensionalAnalysisEquations
    from .equations.cmb_observables import CmbObservablesEquations
except ImportError:
    # Fallback for direct execution
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'equations'))
    from schwarzschild import schwarzschildEquations
    from flux_law import FluxLawEquations
    from adm import AdmEquations
    from kerr_equations import KerrEquations
    from cosmological import CosmologicalEquations
    from christoffel import ChristoffelEquations
    from conservation import ConservationEquations
    from coordinate_transformations import CoordinateTransformationEquations
    from reissner_nordstrom import ReissnerNordstromEquations
    from vaidya import VaidyaEquations
    from kerr_newman import KerrNewmanEquations
    from dimensional_analysis import DimensionalAnalysisEquations
    from cmb_observables import CmbObservablesEquations

class SymbolicGATGEquations:
    """
    Symbolic equation repository that connects documentation with computation

    Now organized with category-specific modules for better maintainability.
    Each equation entry contains:
    - symbolic: The actual SageMath symbolic expression
    - string: Human-readable string representation
    - latex: LaTeX formatted version
    - description: Physical/mathematical significance
    """

    def __init__(self):
        """Initialize symbolic variables and equation definitions"""
        self.vars = {}  # Variable registry
        self.equations = {}  # Equation storage
        self.categories = {
            'schwarzschild': [],
            'kerr': [],
            'flux_law': [],
            'adm': [],
            'cosmological': [],
            'christoffel': [],
            'conservation': [],
            'coordinate_transformations': [],
            'reissner_nordstrom': [],
            'vaidya': [],
            'kerr_newman': [],
            'dimensional_analysis': [],
            'cmb_observables': []
        }
        self.category_modules = {}

        # Setup all variables and equations
        self._setup_variables()
        self._load_category_modules()

    def _setup_variables(self):
        """Setup consistent symbolic variables used across all equations"""
        # Coordinates
        self.vars['t'] = var('t', domain='real')
        self.vars['r'] = var('r', domain='positive')
        self.vars['theta'] = var('theta', domain='real')
        self.vars['phi'] = var('phi', domain='real')

        # Shorter aliases
        self.vars['th'] = self.vars['theta']
        self.vars['ph'] = self.vars['phi']

        # Physical parameters
        self.vars['M'] = var('M', domain='positive')  # Mass
        self.vars['G'] = var('G', domain='positive')  # Gravitational constant
        self.vars['c'] = var('c', domain='positive')  # Speed of light
        self.vars['a'] = var('a', domain='real')  # Kerr parameter
        self.vars['Q'] = var('Q', domain='real')  # Electric charge
        self.vars['Lambda'] = var('Lambda', domain='real')  # Cosmological constant

        # Derived quantities
        self.vars['r_s'] = var('r_s', domain='positive')  # Schwarzschild radius
        self.vars['rs'] = self.vars['r_s']  # Alias

        # Cosmological variables
        self.vars['H'] = var('H', domain='real')  # Hubble parameter
        self.vars['H_0'] = var('H_0', domain='positive')  # Hubble constant
        self.vars['k'] = var('k', domain='real')  # Curvature parameter (-1, 0, +1)
        self.vars['z'] = var('z', domain='positive')  # Redshift
        self.vars['rho'] = var('rho', domain='positive')  # Energy density
        self.vars['p'] = var('p', domain='real')  # Pressure
        self.vars['rho_m'] = var('rho_m', domain='positive')  # Matter density
        self.vars['rho_r'] = var('rho_r', domain='positive')  # Radiation density
        self.vars['rho_Lambda'] = var('rho_Lambda', domain='positive')  # Dark energy density
        self.vars['Omega_m'] = var('Omega_m', domain='positive')  # Matter density parameter
        self.vars['Omega_Lambda'] = var('Omega_Lambda', domain='positive')  # Dark energy parameter
        self.vars['w'] = var('w', domain='real')  # Equation of state parameter

        # Apply physical assumptions
        assume(self.vars['theta'] >= 0)
        assume(self.vars['theta'] <= pi)
        assume(self.vars['phi'] >= 0)
        assume(self.vars['phi'] <= 2*pi)
        assume(self.vars['k'] >= -1)
        assume(self.vars['k'] <= 1)

    def _load_category_modules(self):
        """Load all category modules and their equations"""

        # Initialize category modules
        self.category_modules['schwarzschild'] = schwarzschildEquations(self.vars)
        self.category_modules['flux_law'] = FluxLawEquations(self.vars)
        self.category_modules['adm'] = AdmEquations(self.vars)
        self.category_modules['kerr'] = KerrEquations(self.vars)
        self.category_modules['cosmological'] = CosmologicalEquations(self.vars)
        self.category_modules['christoffel'] = ChristoffelEquations(self.vars)
        self.category_modules['conservation'] = ConservationEquations(self.vars)
        self.category_modules['coordinate_transformations'] = CoordinateTransformationEquations(self.vars)
        self.category_modules['reissner_nordstrom'] = ReissnerNordstromEquations(self.vars)
        self.category_modules['vaidya'] = VaidyaEquations(self.vars)
        self.category_modules['kerr_newman'] = KerrNewmanEquations(self.vars)
        self.category_modules['dimensional_analysis'] = DimensionalAnalysisEquations(self.vars)
        self.category_modules['cmb_observables'] = CmbObservablesEquations(self.vars)

        # Load equations from each module
        for category_name, module in self.category_modules.items():
            category_equations = module.get_equations()
            self.equations.update(category_equations)
            self.categories[category_name] = module.get_equation_keys()

    # Access methods
    def get_symbolic(self, equation_key: str):
        """Return the symbolic expression for an equation"""
        if equation_key not in self.equations:
            raise KeyError(f"Equation '{equation_key}' not found")
        return self.equations[equation_key]['symbolic']

    def get_string(self, equation_key: str):
        """Return the string representation for an equation"""
        if equation_key not in self.equations:
            raise KeyError(f"Equation '{equation_key}' not found")
        return self.equations[equation_key]['string']

    def get_latex(self, equation_key: str):
        """Return the LaTeX representation for an equation"""
        if equation_key not in self.equations:
            raise KeyError(f"Equation '{equation_key}' not found")
        return self.equations[equation_key]['latex']

    def get_description(self, equation_key: str):
        """Return the description for an equation"""
        if equation_key not in self.equations:
            raise KeyError(f"Equation '{equation_key}' not found")
        return self.equations[equation_key]['description']

    # Verification methods - delegate to category modules
    def verify_schwarzschild_ode(self) -> bool:
        return self.category_modules['schwarzschild'].verify_schwarzschild_ode()

    def verify_flux_law_consistency(self) -> bool:
        return self.category_modules['flux_law'].verify_flux_law_consistency()

    def verify_kerr_schwarzschild_limit(self) -> bool:
        return self.category_modules['kerr'].verify_kerr_schwarzschild_limit()

    def verify_kerr_horizons(self) -> bool:
        return self.category_modules['kerr'].verify_kerr_horizons()

    def verify_kerr_metric_consistency(self) -> bool:
        return self.category_modules['kerr'].verify_kerr_metric_consistency()

    def verify_lapse_first_hubble_relation(self) -> bool:
        return self.category_modules['cosmological'].verify_lapse_first_hubble_relation()

    def verify_matter_scaling(self) -> bool:
        return self.category_modules['cosmological'].verify_matter_scaling()

    def verify_radiation_scaling(self) -> bool:
        return self.category_modules['cosmological'].verify_radiation_scaling()

    def verify_flatness_condition(self) -> bool:
        return self.category_modules['cosmological'].verify_flatness_condition()

    def verify_cosmological_redshift_relation(self) -> bool:
        return self.category_modules['cosmological'].verify_cosmological_redshift_relation()

    def verify_christoffel_symmetry(self) -> bool:
        return self.category_modules['christoffel'].verify_christoffel_symmetry()

    def verify_christoffel_connection_property(self) -> bool:
        return self.category_modules['christoffel'].verify_christoffel_connection_property()

    def verify_christoffel_spherical_symmetry(self) -> bool:
        return self.category_modules['christoffel'].verify_christoffel_spherical_symmetry()

    def verify_stress_energy_conservation_structure(self) -> bool:
        return self.category_modules['conservation'].verify_stress_energy_conservation_structure()

    def verify_covariant_divergence_consistency(self) -> bool:
        return self.category_modules['conservation'].verify_covariant_divergence_consistency()

    def verify_momentum_constraint_cancellation(self) -> bool:
        return self.category_modules['conservation'].verify_momentum_constraint_cancellation()

    def verify_coordinate_transformation_consistency(self) -> bool:
        return self.category_modules['coordinate_transformations'].verify_coordinate_transformation_consistency()

    def verify_vaidya_transformation_structure(self) -> bool:
        return self.category_modules['coordinate_transformations'].verify_vaidya_transformation_structure()

    def verify_reissner_nordstrom_horizons(self) -> bool:
        return self.category_modules['reissner_nordstrom'].verify_reissner_nordstrom_horizons()

    def verify_reissner_nordstrom_schwarzschild_limit(self) -> bool:
        return self.category_modules['reissner_nordstrom'].verify_schwarzschild_limit()

    def verify_vaidya_limiting_cases(self) -> bool:
        return self.category_modules['vaidya'].verify_vaidya_limiting_cases()

    def verify_vaidya_mass_function_consistency(self) -> bool:
        return self.category_modules['vaidya'].verify_mass_function_consistency()

    def verify_kerr_newman_horizons(self) -> bool:
        return self.category_modules['kerr_newman'].verify_kerr_newman_horizons()

    def verify_kerr_newman_limits(self) -> bool:
        return self.category_modules['kerr_newman'].verify_kerr_newman_limits()

    def verify_flux_law_dimensional_consistency(self) -> bool:
        return self.category_modules['dimensional_analysis'].verify_flux_law_dimensional_consistency()

    def verify_cmb_recombination_physics(self) -> bool:
        return self.category_modules['cmb_observables'].verify_recombination_physics()

    def evaluate(self, equation_key: str, substitutions: Dict[Any, Any]):
        """Evaluate an equation with given substitutions"""
        expr = self.get_symbolic(equation_key)
        return expr.substitute(substitutions).simplify_full()

    def list_equations(self, category: Optional[str] = None):
        """List all equations or equations in a specific category"""
        if category:
            if category not in self.categories:
                raise KeyError(f"Category '{category}' not found")
            return self.categories[category]
        return list(self.equations.keys())

    def get_equation_info(self, equation_key: str) -> Dict[str, Any]:
        """Return complete information about an equation"""
        if equation_key not in self.equations:
            raise KeyError(f"Equation '{equation_key}' not found")
        return self.equations[equation_key].copy()

    def generate_latex_document(self, category: Optional[str] = None) -> str:
        """Generate a LaTeX document fragment for equations"""
        equations_to_process = (
            self.categories[category] if category
            else list(self.equations.keys())
        )

        latex_lines = []
        for eq_key in equations_to_process:
            eq = self.equations[eq_key]
            latex_lines.append(f"% {eq['description']}")
            latex_lines.append(f"\\begin{{equation}}")
            latex_lines.append(f"  {eq['latex']}")
            latex_lines.append(f"  \\label{{eq:{eq_key}}}")
            latex_lines.append(f"\\end{{equation}}")
            latex_lines.append("")

        return "\n".join(latex_lines)

# Create module-level instance for easy import
equations = SymbolicGATGEquations()

# Export convenience functions
def get_symbolic(equation_key: str):
    """Get symbolic expression for an equation"""
    return equations.get_symbolic(equation_key)

def get_string(equation_key: str):
    """Get string representation for an equation"""
    return equations.get_string(equation_key)

def verify_schwarzschild():
    """Verify Schwarzschild equations"""
    return equations.verify_schwarzschild_ode()

def verify_flux_law():
    """Verify flux law consistency"""
    return equations.verify_flux_law_consistency()

def verify_kerr_schwarzschild_limit():
    """Verify Kerr reduces to Schwarzschild when a → 0"""
    return equations.verify_kerr_schwarzschild_limit()

def verify_kerr_horizons():
    """Verify Kerr event horizons satisfy Δ = 0"""
    return equations.verify_kerr_horizons()

def verify_kerr_consistency():
    """Verify internal consistency of Kerr auxiliary functions"""
    return equations.verify_kerr_metric_consistency()

def verify_lapse_first_hubble():
    """Verify consistency between standard and lapse-first Hubble relations"""
    return equations.verify_lapse_first_hubble_relation()

def verify_matter_scaling():
    """Verify matter density scaling ρ ∝ a^(-3)"""
    return equations.verify_matter_scaling()

def verify_radiation_scaling():
    """Verify radiation density scaling ρ ∝ a^(-4)"""
    return equations.verify_radiation_scaling()

def verify_flatness_condition():
    """Verify density parameters sum to unity for flat universe"""
    return equations.verify_flatness_condition()

def verify_cosmological_redshift():
    """Verify redshift-scale factor relation"""
    return equations.verify_cosmological_redshift_relation()

def verify_christoffel_symmetry():
    """Verify Christoffel symbols satisfy expected symmetries"""
    return equations.verify_christoffel_symmetry()

def verify_christoffel_connection():
    """Verify Christoffel symbols satisfy the connection property"""
    return equations.verify_christoffel_connection_property()

def verify_christoffel_spherical_symmetry():
    """Verify Christoffel symbols respect spherical symmetry"""
    return equations.verify_christoffel_spherical_symmetry()

# Export main class and functions
__all__ = [
    'SymbolicGATGEquations',
    'symbolic_equations',
    'get_symbolic',
    'get_string',
    'verify_schwarzschild',
    'verify_flux_law',
    'verify_kerr_schwarzschild_limit',
    'verify_kerr_horizons',
    'verify_kerr_consistency',
    'verify_lapse_first_hubble',
    'verify_matter_scaling',
    'verify_radiation_scaling',
    'verify_flatness_condition',
    'verify_cosmological_redshift',
    'verify_christoffel_symmetry',
    'verify_christoffel_connection',
    'verify_christoffel_spherical_symmetry'
]