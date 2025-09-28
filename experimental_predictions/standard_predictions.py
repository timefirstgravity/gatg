#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Standard GR Experimental Predictions

Implements standard GR approach to experimental predictions and observable consequences
Computes actual symbolic expressions for gravitational effects and measurable quantities.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

def compute_perihelion_precession():
    """
    Compute perihelion precession from standard GR

    Returns:
        dict: Perihelion precession calculation
    """
    # Orbital parameters
    a = var('a')  # Semi-major axis
    e = var('e')  # Eccentricity
    c = var('c')  # Speed of light
    G = var('G')  # Gravitational constant
    M = var('M')  # Central mass

    # Standard GR prediction for perihelion advance per orbit
    # Δφ = 6πGM/(ac²(1-e²))
    delta_phi_per_orbit = 6*pi*G*M/(a*c**2*(1-e**2))

    # Convert to arcseconds per century for comparison with observations
    # Using orbital period T = 2π√(a³/GM)
    T = 2*pi*sqrt(a**3/(G*M))
    orbits_per_century = var('century')/T  # century = 100 years in seconds

    precession_per_century = delta_phi_per_orbit * orbits_per_century

    return {
        'perihelion_advance_per_orbit': delta_phi_per_orbit,
        'orbital_period': T,
        'precession_per_century': precession_per_century,
        'formula': 'Δφ = 6πGM/(ac²(1-e²))',
        'physical_origin': 'Spacetime curvature around massive body'
    }

def compute_light_deflection():
    """
    Compute light deflection from standard GR

    Returns:
        dict: Light deflection calculation
    """
    # Impact parameter and source mass
    b = var('b')  # Impact parameter (closest approach distance)
    M = var('M')  # Deflecting mass
    G = var('G')
    c = var('c')

    # Standard GR prediction for light deflection angle
    # α = 4GM/(bc²)
    deflection_angle = 4*G*M/(b*c**2)

    # For the Sun: deflection at limb
    R_sun = var('R_sun')  # Solar radius
    M_sun = var('M_sun')  # Solar mass
    deflection_at_limb = 4*G*M_sun/(R_sun*c**2)

    # Express in arcseconds
    # 1 radian ≈ 206265 arcseconds
    arcsec_per_radian = 206265
    deflection_arcsec = deflection_at_limb * arcsec_per_radian

    return {
        'deflection_angle': deflection_angle,
        'deflection_at_solar_limb': deflection_at_limb,
        'deflection_arcseconds': deflection_arcsec,
        'formula': 'α = 4GM/(bc²)',
        'physical_origin': 'Null geodesics in curved spacetime'
    }

def compute_gravitational_redshift():
    """
    Compute gravitational redshift from standard GR

    Returns:
        dict: Gravitational redshift calculation
    """
    # Source and observer positions
    r_source = var('r_s')      # Radial position of source
    r_observer = var('r_o')    # Radial position of observer
    M = var('M')               # Central mass
    G = var('G')
    c = var('c')

    # Schwarzschild metric time dilation factors
    # g_tt(r) = -(1 - 2GM/(rc²))
    lapse_source = sqrt(1 - 2*G*M/(r_source*c**2))
    lapse_observer = sqrt(1 - 2*G*M/(r_observer*c**2))

    # Redshift formula: ν_o/ν_s = lapse_observer/lapse_source
    redshift_factor = lapse_observer/lapse_source

    # For weak field (r >> 2GM/c²): z ≈ GM/(rc²) * (1/r_s - 1/r_o)
    weak_field_redshift = G*M/c**2 * (1/r_source - 1/r_observer)

    return {
        'redshift_factor': redshift_factor,
        'weak_field_approximation': weak_field_redshift,
        'lapse_at_source': lapse_source,
        'lapse_at_observer': lapse_observer,
        'formula': 'ν_o/ν_s = √[g_tt(r_o)/g_tt(r_s)]',
        'physical_origin': 'Time dilation in gravitational field'
    }

def compute_gravitational_waves():
    """
    Compute gravitational wave predictions from standard GR

    Returns:
        dict: Gravitational wave properties
    """
    # Wave parameters
    h = var('h')           # Wave amplitude
    f = var('f')           # Frequency
    c = var('c')
    G = var('G')

    # Source parameters
    M1 = var('M1')         # Mass 1
    M2 = var('M2')         # Mass 2
    r = var('r')           # Distance to source

    # Chirp mass for binary system
    M_chirp = (M1*M2)**(3/5) / (M1 + M2)**(1/5)

    # Characteristic strain amplitude
    # h ~ (G/c⁴) * (M_chirp*f²) / r
    strain_amplitude = (G/c**4) * M_chirp * f**2 / r

    # Frequency evolution (chirp)
    # df/dt = (96π/5) * (πGM_chirp/c³)^(5/3) * f^(11/3)
    frequency_evolution = (96*pi/5) * (pi*G*M_chirp/c**3)**(5/3) * f**(11/3)

    # Energy loss rate
    # dE/dt = -(32/5) * G⁴/(c⁵) * (M1*M2)²*(M1+M2)/r⁵
    energy_loss_rate = -(32/5) * G**4/c**5 * (M1*M2)**2 * (M1+M2) / r**5

    return {
        'strain_amplitude': strain_amplitude,
        'frequency_evolution': frequency_evolution,
        'energy_loss_rate': energy_loss_rate,
        'chirp_mass': M_chirp,
        'wave_speed': c,
        'physical_origin': 'Spacetime ripples from accelerating masses'
    }

def compute_time_delay():
    """
    Compute Shapiro time delay from standard GR

    Returns:
        dict: Time delay calculation
    """
    # Signal path parameters
    r_source = var('r_s')      # Source distance
    r_observer = var('r_o')    # Observer distance
    r_min = var('r_min')       # Minimum approach distance
    M = var('M')               # Deflecting mass
    G = var('G')
    c = var('c')

    # Time delay formula for signal passing near massive body
    # Δt = (4GM/c³) * ln[(r_s + r_o + d)/(r_s + r_o - d)]
    # where d is the distance between source and observer
    d = sqrt((r_source - r_observer)**2)  # Simplified

    # Shapiro delay approximation
    time_delay = (4*G*M/c**3) * ln((r_source + r_observer + d)/(r_source + r_observer - d))

    # For small deflection: Δt ≈ (2GM/c³) * ln(4*r_s*r_o/r_min²)
    small_angle_delay = (2*G*M/c**3) * ln(4*r_source*r_observer/r_min**2)

    return {
        'time_delay': time_delay,
        'small_angle_approximation': small_angle_delay,
        'formula': 'Δt = (4GM/c³) ln[(r_s+r_o+d)/(r_s+r_o-d)]',
        'physical_origin': 'Light follows curved spacetime geodesics'
    }