#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
ADM 3+1 Decomposition Functions for GATG
Fundamental ADM formalism operations that are reusable across modules
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

def compute_extrinsic_curvature(gamma_rr, gamma_thth, gamma_phph, N, t):
    """Compute extrinsic curvature tensor K_ij"""
    # Time derivatives
    dt_gamma_rr = diff(gamma_rr, t)
    dt_gamma_thth = diff(gamma_thth, t)
    dt_gamma_phph = diff(gamma_phph, t)

    # K_ij components
    K_rr = -(1/(2*N)) * dt_gamma_rr
    K_thth = -(1/(2*N)) * dt_gamma_thth
    K_phph = -(1/(2*N)) * dt_gamma_phph

    return K_rr, K_thth, K_phph, dt_gamma_rr, dt_gamma_thth, dt_gamma_phph

def compute_mixed_extrinsic_curvature(K_rr, K_thth, K_phph, A, r, th):
    """Compute mixed extrinsic curvature K^i_j"""
    # K^i_j = γ^ik K_kj where γ^ik is the inverse 3-metric
    # Inverse 3-metric components:
    # γ^rr = A = e^(2Φ)
    # γ^θθ = 1/r²
    # γ^φφ = 1/(r² sin²θ)

    Krr_up = A * K_rr  # K^r_r = γ^rr K_rr
    Kthth_up = (1/r**2) * K_thth  # K^θ_θ = γ^θθ K_θθ
    Kphph_up = (1/(r**2 * sin(th)**2)) * K_phph  # K^φ_φ = γ^φφ K_φφ

    # Trace K = K^i_i = K^r_r + K^θ_θ + K^φ_φ
    K_trace = Krr_up + Kthth_up + Kphph_up
    return Krr_up, Kthth_up, Kphph_up, K_trace

def compute_Y_tensor(Krr_up, K_trace):
    """Compute Y tensor Y^j_i = K^j_i - δ^j_i K"""
    Y_rr = Krr_up - K_trace  # Should be 0
    Y_thth = -K_trace  # Y^θ_θ = K^θ_θ - K = 0 - K = -K
    Y_phph = -K_trace  # Y^φ_φ = K^φ_φ - K = 0 - K = -K
    return Y_rr, Y_thth, Y_phph

def compute_covariant_divergence(Y_thth, Y_phph, r):
    """Compute covariant divergence D_j Y^j_r"""
    # Christoffel symbols for spherical metric
    Gamma_theta_rtheta = 1/r
    Gamma_phi_rphi = 1/r

    # D_j Y^j_r = -Γ^m_rj Y^j_m (since Y^j_r = 0 for all j)
    DY = -Gamma_theta_rtheta * Y_thth - Gamma_phi_rphi * Y_phph
    return DY, Gamma_theta_rtheta, Gamma_phi_rphi