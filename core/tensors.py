#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Core: Tensor Computations
Common tensor operations for all GATG derivations
"""

from sage.all import *

def compute_ricci_tensor(g):
    """
    Compute Ricci tensor and scalar for given metric
    """
    R = g.ricci()
    Rs = g.ricci_scalar()
    return R, Rs

def compute_einstein_tensor(g):
    """
    Compute Einstein tensor G_μν = R_μν - (1/2) R g_μν
    """
    R = g.ricci()
    Rs = g.ricci_scalar()
    G = R - (Rs/2)*g
    return G, R, Rs

def extract_einstein_components(G, simplify=True, include_off_diagonal=False):
    """
    Extract and optionally simplify Einstein tensor components

    Parameters:
        G: Einstein tensor
        simplify: Whether to apply simplify_full() to components
        include_off_diagonal: Whether to include off-diagonal components (needed for Kerr, etc.)

    Returns:
        Dictionary of components with appropriate labels
    """
    components = {}

    # Always include diagonal components
    components['tt'] = G[0,0].expr()
    components['rr'] = G[1,1].expr()
    components['thth'] = G[2,2].expr()
    components['phph'] = G[3,3].expr()

    # Include off-diagonal components if requested (for non-spherically symmetric cases)
    if include_off_diagonal:
        components['tr'] = G[0,1].expr()
        components['tth'] = G[0,2].expr()
        components['tph'] = G[0,3].expr()
        components['rth'] = G[1,2].expr()
        components['rph'] = G[1,3].expr()
        components['thph'] = G[2,3].expr()

    if simplify:
        for key in components:
            components[key] = components[key].simplify_full()

    return components


def verify_vacuum_solution(G, tolerance=0, include_off_diagonal=False):
    """
    Verify that Einstein tensor vanishes (vacuum solution)

    Parameters:
        G: Einstein tensor
        tolerance: Numerical tolerance (currently unused)
        include_off_diagonal: Whether to check off-diagonal components too

    Returns:
        Tuple: (all_zero, components)
    """
    components = extract_einstein_components(G, simplify=True, include_off_diagonal=include_off_diagonal)

    all_zero = True
    for key, value in components.items():
        if hasattr(value, 'simplify_full'):
            simplified = value.simplify_full()
        else:
            simplified = value

        if simplified != 0:
            all_zero = False

    return all_zero, components

def compute_einstein_tensor_component(M, A, r, t):
    """Compute G^t_r component independently for verification"""
    from .metrics import construct_metric

    g_test = construct_metric(M, A, r, var('theta'))

    # Compute Einstein tensor
    R_test = g_test.ricci()
    Rs_test = g_test.ricci_scalar()
    G_test = R_test - (Rs_test/2)*g_test

    # Get mixed component G^t_r
    G_mixed = G_test.up(g_test)
    Gtr = G_mixed[0,1].expr()

    return Gtr