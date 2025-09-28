#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Core: Metric Construction
Common metric operations for all GATG derivations
"""

from sage.all import *

def construct_metric(M, A, r, th):
    """
    Construct the spherical metric tensor
    ds^2 = -A dt^2 + A^-1 dr^2 + r^2 dtheta^2 + r^2 sin^2(theta) dphi^2
    """
    g = M.metric('g')
    g[0,0] = -A
    g[1,1] = 1/A
    g[2,2] = r**2
    g[3,3] = r**2 * sin(th)**2
    return g

def define_temporal_potential(M, t, r, static=False):
    """
    Define the temporal potential Phi(t,r) and related quantities
    """
    if static:
        Phi = M.scalar_field(function('Phi')(r), name='Phi')
    else:
        Phi = M.scalar_field(function('Phi')(t,r), name='Phi')

    A = exp(2*Phi.expr())
    N = exp(Phi.expr())  # Lapse function
    return Phi, A, N

def define_static_temporal_potential(M, r):
    """
    Define static temporal potential Phi(r) and related quantities
    """
    return define_temporal_potential(M, None, r, static=True)

def compute_3metric_components(A, r, th):
    """
    Compute spatial 3-metric components
    """
    gamma_rr = 1/A
    gamma_thth = r**2
    gamma_phph = r**2 * sin(th)**2
    return gamma_rr, gamma_thth, gamma_phph