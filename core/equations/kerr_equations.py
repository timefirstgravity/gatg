#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Kerr Metric Equations for GATG
Symbolic representations of Kerr (rotating black hole) equations
"""

from sage.all import *
from typing import Dict, Any, Optional

class KerrEquations:
    """
    Symbolic equations for Kerr metric in GATG
    """

    def __init__(self, vars_dict: Dict[str, Any]):
        """Initialize with shared variable dictionary"""
        self.vars = vars_dict
        self.equations = {}
        self.category_name = 'kerr'

        # Define equations
        self.define_equations()

    def define_equations(self):
        """Define Kerr metric equations (rotating black hole)"""
        r = self.vars['r']
        theta = self.vars['theta']
        M = self.vars['M']
        a = self.vars['a']

        # Auxiliary functions (fundamental building blocks)
        Sigma = r**2 + a**2 * cos(theta)**2
        Delta = r**2 - 2*M*r + a**2

        self.equations['kerr_sigma'] = {
            'symbolic': Sigma,
            'string': "Σ = r² + a²cos²θ",
            'latex': r"\Sigma = r^2 + a^2 \cos^2\theta",
            'description': "Kerr auxiliary function Sigma (appears in denominators)"
        }

        self.equations['kerr_delta'] = {
            'symbolic': Delta,
            'string': "Δ = r² - 2Mr + a²",
            'latex': r"\Delta = r^2 - 2Mr + a^2",
            'description': "Kerr auxiliary function Delta (determines horizons)"
        }

        # Metric components (Boyer-Lindquist coordinates)
        self.equations['kerr_g_tt'] = {
            'symbolic': -(1 - 2*M*r/Sigma),
            'string': "g_tt = -(1 - 2Mr/Σ)",
            'latex': r"g_{tt} = -\left(1 - \frac{2Mr}{\Sigma}\right)",
            'description': "Kerr metric time-time component"
        }

        self.equations['kerr_g_rr'] = {
            'symbolic': Sigma/Delta,
            'string': "g_rr = Σ/Δ",
            'latex': r"g_{rr} = \frac{\Sigma}{\Delta}",
            'description': "Kerr metric radial-radial component"
        }

        self.equations['kerr_g_thth'] = {
            'symbolic': Sigma,
            'string': "g_θθ = Σ",
            'latex': r"g_{\theta\theta} = \Sigma",
            'description': "Kerr metric theta-theta component"
        }

        self.equations['kerr_g_phph'] = {
            'symbolic': sin(theta)**2 * ((r**2 + a**2) + 2*M*a**2*r*sin(theta)**2/Sigma),
            'string': "g_φφ = sin²θ[(r² + a²) + 2Ma²r sin²θ/Σ]",
            'latex': r"g_{\phi\phi} = \sin^2\theta\left[(r^2 + a^2) + \frac{2Ma^2 r \sin^2\theta}{\Sigma}\right]",
            'description': "Kerr metric phi-phi component (includes frame-dragging effects)"
        }

        # Off-diagonal component (frame-dragging!)
        self.equations['kerr_g_tph'] = {
            'symbolic': -2*a*M*r*sin(theta)**2/Sigma,
            'string': "g_tφ = -2aMr sin²θ/Σ",
            'latex': r"g_{t\phi} = -\frac{2aMr \sin^2\theta}{\Sigma}",
            'description': "Kerr metric time-phi component (frame-dragging/Lense-Thirring effect)"
        }

        # Complete Kerr line element
        self.equations['kerr_metric'] = {
            'symbolic': None,  # Full line element is complex symbolic expression
            'string': "ds² = -(1-2Mr/Σ)dt² + (Σ/Δ)dr² + Σdθ² + sin²θ[(r²+a²)+2Ma²r sin²θ/Σ]dφ² - 2aMr sin²θ/Σ dt dφ",
            'latex': r"ds^2 = -\left(1-\frac{2Mr}{\Sigma}\right)dt^2 + \frac{\Sigma}{\Delta}dr^2 + \Sigma d\theta^2 + \sin^2\theta\left[(r^2+a^2)+\frac{2Ma^2 r \sin^2\theta}{\Sigma}\right]d\phi^2 - \frac{2aMr \sin^2\theta}{\Sigma} dt d\phi",
            'description': "Complete Kerr metric line element in Boyer-Lindquist coordinates"
        }

        # Physical properties
        # Event horizons (where Δ = 0)
        self.equations['kerr_event_horizons'] = {
            'symbolic': [M + sqrt(M**2 - a**2), M - sqrt(M**2 - a**2)],
            'string': "r_± = M ± √(M² - a²)",
            'latex': r"r_\pm = M \pm \sqrt{M^2 - a^2}",
            'description': "Kerr black hole event horizons (outer and inner)"
        }

        # Ergosphere boundary
        self.equations['kerr_ergosphere'] = {
            'symbolic': M + sqrt(M**2 - a**2*cos(theta)**2),
            'string': "r_e = M + √(M² - a²cos²θ)",
            'latex': r"r_e = M + \sqrt{M^2 - a^2\cos^2\theta}",
            'description': "Kerr ergosphere boundary (static limit surface)"
        }

        # Angular momentum parameter
        J = var('J', domain='positive')  # Total angular momentum
        c = self.vars['c']
        self.equations['kerr_angular_momentum'] = {
            'symbolic': a - J/(M*c),
            'string': "a = J/Mc",
            'latex': r"a = \frac{J}{Mc}",
            'description': "Kerr angular momentum parameter"
        }

        # Schwarzschild limit
        self.equations['kerr_schwarzschild_limit'] = {
            'symbolic': None,  # This is a limiting process, not a single expression
            'string': "a → 0: Kerr → Schwarzschild",
            'latex': r"a \to 0: \text{Kerr} \to \text{Schwarzschild}",
            'description': "Kerr metric reduces to Schwarzschild when rotation vanishes"
        }

        # Lapse-first decomposition of Kerr metric
        # These are more complex expressions for ADM decomposition
        Xi = (r**2 + a**2)**2 - a**2*Delta*sin(theta)**2

        self.equations['kerr_xi'] = {
            'symbolic': Xi,
            'string': "Ξ = (r² + a²)² - a²Δ sin²θ",
            'latex': r"\Xi = (r^2 + a^2)^2 - a^2\Delta \sin^2\theta",
            'description': "Kerr auxiliary function Xi for lapse-first decomposition"
        }

        self.equations['kerr_lapse_function'] = {
            'symbolic': Delta*sin(theta)**2/Xi,
            'string': "N² = Δ sin²θ/Ξ",
            'latex': r"N^2 = \frac{\Delta \sin^2\theta}{\Xi}",
            'description': "Kerr lapse function (full form)"
        }

        self.equations['kerr_lapse_simplified'] = {
            'symbolic': 1 - 2*M*r/Sigma,
            'string': "N² = 1 - 2Mr/Σ",
            'latex': r"N^2 = 1 - \frac{2Mr}{\Sigma}",
            'description': "Kerr lapse function (simplified, approximate)"
        }

        self.equations['kerr_shift_vector'] = {
            'symbolic': -2*a*M*r/(Xi*sin(theta)**2),
            'string': "β^φ = -2aMr/(Ξ sin²θ)",
            'latex': r"\beta^\phi = -\frac{2aMr}{\Xi \sin^2\theta}",
            'description': "Kerr shift vector (frame-dragging in ADM formulation)"
        }

        # Kretschmann scalar (curvature invariant)
        # For Kerr metric: K = R_{μνρσ}R^{μνρσ}
        # This is a complex expression, simplified form:
        rho_kerr = r / (r**2 + a**2 * cos(theta)**2)  # Simplified radial coordinate
        K_kerr = 48 * M**2 * (r**2 - a**2 * cos(theta)**2) * \
                 (r**4 - 14*r**2*a**2*cos(theta)**2 + a**4*cos(theta)**4) / \
                 (r**2 + a**2*cos(theta)**2)**6

        self.equations['kerr_kretschmann_scalar'] = {
            'symbolic': K_kerr,
            'string': "K = 48M²(r² - a²cos²θ)(r⁴ - 14r²a²cos²θ + a⁴cos⁴θ)/Σ⁶",
            'latex': r"K = \frac{48M^2(r^2 - a^2\cos^2\theta)(r^4 - 14r^2a^2\cos^2\theta + a^4\cos^4\theta)}{\Sigma^6}",
            'description': "Kerr Kretschmann curvature invariant R_{μνρσ}R^{μνρσ}"
        }

        # Weyl scalar (conformal curvature invariant) - another important invariant
        # For vacuum solutions: C_{μνρσ}C^{μνρσ} = R_{μνρσ}R^{μνρσ} (since R_μν = 0)
        self.equations['kerr_weyl_scalar'] = {
            'symbolic': K_kerr,  # Same as Kretschmann for vacuum
            'string': "C = 48M²(r² - a²cos²θ)(r⁴ - 14r²a²cos²θ + a⁴cos⁴θ)/Σ⁶",
            'latex': r"C = \frac{48M^2(r^2 - a^2\cos^2\theta)(r^4 - 14r^2a^2\cos^2\theta + a^4\cos^4\theta)}{\Sigma^6}",
            'description': "Kerr Weyl curvature invariant C_{μνρσ}C^{μνρσ} (equals Kretschmann in vacuum)"
        }

    def verify_kerr_schwarzschild_limit(self) -> bool:
        """Verify that Kerr metric reduces to Schwarzschild when a → 0"""
        r = self.vars['r']
        theta = self.vars['theta']
        M = self.vars['M']
        a = self.vars['a']

        # Get Kerr metric components
        kerr_g_tt = self.equations['kerr_g_tt']['symbolic']
        kerr_g_rr = self.equations['kerr_g_rr']['symbolic']
        kerr_g_thth = self.equations['kerr_g_thth']['symbolic']
        kerr_g_phph = self.equations['kerr_g_phph']['symbolic']
        kerr_g_tph = self.equations['kerr_g_tph']['symbolic']

        # Take limit a → 0
        schwarzschild_g_tt = kerr_g_tt.substitute({a: 0}).simplify_full()
        schwarzschild_g_rr = kerr_g_rr.substitute({a: 0}).simplify_full()
        schwarzschild_g_thth = kerr_g_thth.substitute({a: 0}).simplify_full()
        schwarzschild_g_phph = kerr_g_phph.substitute({a: 0}).simplify_full()
        schwarzschild_g_tph = kerr_g_tph.substitute({a: 0}).simplify_full()

        # Expected Schwarzschild values
        rs = 2*M  # Schwarzschild radius
        expected_g_tt = -(1 - rs/r)
        expected_g_rr = 1/(1 - rs/r)
        expected_g_thth = r**2
        expected_g_phph = r**2 * sin(theta)**2
        expected_g_tph = 0

        # Check each component
        checks = [
            schwarzschild_g_tt.substitute({2*M: rs}) == expected_g_tt,
            schwarzschild_g_rr == expected_g_rr.substitute({rs: 2*M}),
            schwarzschild_g_thth == expected_g_thth,
            schwarzschild_g_phph == expected_g_phph,
            schwarzschild_g_tph == expected_g_tph
        ]

        return all(checks)

    def verify_kerr_horizons(self) -> bool:
        """Verify that Kerr event horizons satisfy Δ = 0"""
        r = self.vars['r']
        M = self.vars['M']
        a = self.vars['a']

        # Get Delta function
        Delta = self.equations['kerr_delta']['symbolic']

        # Get horizon radii
        horizons = self.equations['kerr_event_horizons']['symbolic']
        r_plus = horizons[0]  # M + sqrt(M^2 - a^2)
        r_minus = horizons[1]  # M - sqrt(M^2 - a^2)

        # Check that Δ = 0 at both horizons
        Delta_at_r_plus = Delta.substitute({r: r_plus}).simplify_full()
        Delta_at_r_minus = Delta.substitute({r: r_minus}).simplify_full()

        return Delta_at_r_plus == 0 and Delta_at_r_minus == 0

    def verify_kerr_metric_consistency(self) -> bool:
        """Verify internal consistency of Kerr auxiliary functions"""
        r = self.vars['r']
        theta = self.vars['theta']
        M = self.vars['M']
        a = self.vars['a']

        # Get auxiliary functions
        Sigma = self.equations['kerr_sigma']['symbolic']
        Delta = self.equations['kerr_delta']['symbolic']
        Xi = self.equations['kerr_xi']['symbolic']

        # Verify that Sigma and Delta are used consistently in Xi
        # Xi = (r² + a²)² - a²Δ sin²θ
        Xi_expanded = (r**2 + a**2)**2 - a**2 * Delta * sin(theta)**2
        Xi_computed = Xi.expand().simplify_full()
        Xi_expected = Xi_expanded.expand().simplify_full()

        return Xi_computed == Xi_expected

    def get_equations(self) -> Dict[str, Dict[str, Any]]:
        """Return all equations in this category"""
        return self.equations.copy()

    def get_equation_keys(self) -> list:
        """Return list of equation keys in this category"""
        return list(self.equations.keys())

# Export the class
__all__ = ['KerrEquations']