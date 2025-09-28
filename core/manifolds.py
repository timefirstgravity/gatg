#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Core: Manifold and Coordinate Setup
Common manifold operations for all GATG derivations
"""

from sage.all import *

def setup_manifold_and_coordinates(r_range='(0,+oo)'):
    """
    Setup 4D Lorentzian manifold with spherical coordinates

    Args:
        r_range: String specifying r coordinate range
                 '(0,+oo)' for Schwarzschild (default)
                 '(-oo,+oo)' for Boyer-Lindquist (Kerr)

    Returns:
        M: 4D Lorentzian manifold
        X: Spherical coordinate chart (t, r, theta, phi)
        t, r, th, ph: Coordinate functions
    """
    M = Manifold(4, 'M', structure='Lorentzian')
    chart_spec = f't r:{r_range} th:(0,pi):theta ph:(0,2*pi):phi'
    X = M.chart(chart_spec)
    t, r, th, ph = X[:]
    return M, X, t, r, th, ph

def setup_static_manifold_and_coordinates():
    """
    Setup 4D Lorentzian manifold with spherical coordinates for static case
    Same as setup_manifold_and_coordinates but with clearer naming for static scenarios

    Returns:
        M: 4D Lorentzian manifold
        X: Spherical coordinate chart (t, r, theta, phi)
        t, r, th, ph: Coordinate functions
    """
    return setup_manifold_and_coordinates()