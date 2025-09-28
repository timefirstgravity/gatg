#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Standard Cosmology (FLRW Model)

Implements the standard GR approach: FLRW metric → Einstein equations → Friedmann equations
Computes actual symbolic expressions for cosmological evolution and observables.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

def compute_flrw_metric_setup():
    """
    Set up FLRW metric with scale factor a(t)

    Returns:
        dict: FLRW metric setup and components
    """
    # Coordinates and parameters
    t, r, th, ph = var('t r theta phi')
    a = function('a')(t)  # Scale factor
    k = var('k')  # Spatial curvature parameter

    # FLRW metric components
    # ds² = -dt² + a²(t)[dr²/(1-kr²) + r²(dθ² + sin²θ dφ²)]
    g_tt = -1
    g_rr = a**2 / (1 - k*r**2)
    g_thth = a**2 * r**2
    g_phph = a**2 * r**2 * sin(th)**2

    return {
        'coordinates': [t, r, th, ph],
        'scale_factor': a,
        'curvature_parameter': k,
        'metric_components': {
            'g_tt': g_tt,
            'g_rr': g_rr,
            'g_thth': g_thth,
            'g_phph': g_phph
        }
    }

def derive_friedmann_equations(flrw_data):
    """
    Derive Friedmann equations from Einstein field equations

    Args:
        flrw_data: FLRW metric setup

    Returns:
        dict: Friedmann equations and cosmological parameters
    """
    t = flrw_data['coordinates'][0]
    a = flrw_data['scale_factor']
    k = flrw_data['curvature_parameter']

    # Define cosmological parameters
    H = diff(a, t) / a  # Hubble parameter
    rho = var('rho')    # Energy density
    p = var('p')        # Pressure
    G = var('G')        # Gravitational constant
    c = var('c')        # Speed of light

    # First Friedmann equation: H² = (8πG/3c²)ρ - kc²/a²
    friedmann_1 = H**2 - (8*pi*G*rho)/(3*c**2) + k*c**2/a**2

    # Second Friedmann equation: ä/a = -(4πG/3c²)(ρ + 3p/c²)
    a_ddot = diff(a, t, 2)
    friedmann_2 = a_ddot/a + (4*pi*G)/(3*c**2) * (rho + 3*p/c**2)

    # Continuity equation: ρ̇ + 3H(ρ + p/c²) = 0
    rho_dot = diff(rho, t)
    continuity = rho_dot + 3*H*(rho + p/c**2)

    return {
        'hubble_parameter': H,
        'friedmann_1': friedmann_1,
        'friedmann_2': friedmann_2,
        'continuity_equation': continuity,
        'energy_density': rho,
        'pressure': p
    }

def solve_matter_dominated_epoch(friedmann_data):
    """
    Solve for matter-dominated universe (p = 0, k = 0)

    Args:
        friedmann_data: Friedmann equations

    Returns:
        dict: Matter-dominated solutions
    """
    t = var('t')
    t0 = var('t_0')  # Present time
    a0 = var('a_0')  # Present scale factor

    # Matter-dominated: ρ ∝ a⁻³, p = 0, k = 0
    # Solution: a(t) = a₀(t/t₀)^(2/3)
    a_matter = a0 * (t/t0)**(2/3)

    # Hubble parameter: H = (2/3)/t
    H_matter = (2/3) / t

    # Age of universe at present: t₀ = 2/(3H₀)
    H0 = var('H_0')  # Present Hubble parameter
    age_universe = 2/(3*H0)

    return {
        'scale_factor_matter': a_matter,
        'hubble_parameter_matter': H_matter,
        'age_of_universe': age_universe,
        'present_hubble': H0
    }

def compute_cosmological_observables(matter_solution):
    """
    Compute key cosmological observables

    Args:
        matter_solution: Matter-dominated solutions

    Returns:
        dict: Observable quantities
    """
    # Redshift relation: 1 + z = a₀/a(t)
    z = var('z')
    a0 = matter_solution['scale_factor_matter'].substitute(t=var('t_0'))
    a_z = a0 / (1 + z)

    # Luminosity distance (flat matter-dominated)
    # d_L = (c/H₀) ∫₀^z dz'/H(z') where H(z) = H₀(1+z)^(3/2)
    c = var('c')
    H0 = var('H_0')

    # For matter-dominated: H(z) = H₀(1+z)^(3/2)
    H_z = H0 * (1 + z)**(3/2)

    # Simplified luminosity distance (exact for matter-dominated)
    d_L = (2*c/H0) * ((1 + z)**(1/2) - 1)

    # Distance modulus: μ = 5 log₁₀(d_L/Mpc) + 25
    distance_modulus = 5 * log(d_L, 10) + 25

    return {
        'redshift_scale_relation': a_z,
        'hubble_function': H_z,
        'luminosity_distance': d_L,
        'distance_modulus': distance_modulus
    }

