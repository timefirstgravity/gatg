#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
ADM Constraint Computations for GATG
Functions for computing Hamiltonian and momentum constraints in the ADM formalism
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sage.all import *

def compute_general_lambda_momentum_constraint(M, t, r):
    """
    Compute momentum constraint for general metric ds² = -e^(2Φ)dt² + e^(2Λ)dr² + r²dΩ²
    WITHOUT assuming Λ = -Φ, to show the condition under which our flux law emerges
    """
    # Define general temporal potential functions (local variables)
    Phi_general = function('Phi')(t, r)
    Lambda_general = function('Lambda')(t, r)

    # General 3-metric components for ds² = -e^(2Φ)dt² + e^(2Λ)dr² + r²dΩ²
    gamma_rr_general = exp(2*Lambda_general)
    gamma_thth_general = r**2
    gamma_phph_general = r**2 * sin(var('theta'))**2

    # Lapse function
    N_general = exp(Phi_general)

    # Compute time derivatives (general case)
    dt_gamma_rr_general = diff(gamma_rr_general, t)
    dt_gamma_thth_general = diff(gamma_thth_general, t)  # = 0
    dt_gamma_phph_general = diff(gamma_phph_general, t)  # = 0

    # Extrinsic curvature components K_ij = -(1/2N) ∂_t γ_ij
    K_rr_general = -(1/(2*N_general)) * dt_gamma_rr_general
    K_thth_general = -(1/(2*N_general)) * dt_gamma_thth_general  # = 0
    K_phph_general = -(1/(2*N_general)) * dt_gamma_phph_general  # = 0

    # Simplify K_rr
    K_rr_general = K_rr_general.simplify_full()

    # Inverse 3-metric components
    gamma_rr_inv_general = exp(-2*Lambda_general)
    gamma_thth_inv_general = 1/r**2
    gamma_phph_inv_general = 1/(r**2 * sin(var('theta'))**2)

    # Mixed extrinsic curvature K^i_j = γ^ik K_kj
    Krr_up_general = gamma_rr_inv_general * K_rr_general
    Kthth_up_general = gamma_thth_inv_general * K_thth_general  # = 0
    Kphph_up_general = gamma_phph_inv_general * K_phph_general  # = 0

    # Trace K = K^i_i
    K_trace_general = Krr_up_general + Kthth_up_general + Kphph_up_general
    K_trace_general = K_trace_general.simplify_full()

    # Y tensor Y^j_i = K^j_i - δ^j_i K
    Y_rr_general = Krr_up_general - K_trace_general  # Should be 0
    Y_thth_general = Kthth_up_general - K_trace_general  # = -K_trace
    Y_phph_general = Kphph_up_general - K_trace_general  # = -K_trace

    # Covariant divergence D_j Y^j_r for general case
    # Only angular terms contribute: -Γ^θ_rθ Y^θ_θ - Γ^φ_rφ Y^φ_φ
    Gamma_theta_rtheta_general = 1/r  # Same as before
    Gamma_phi_rphi_general = 1/r      # Same as before

    DY_general = -Gamma_theta_rtheta_general * Y_thth_general - Gamma_phi_rphi_general * Y_phph_general
    DY_general = DY_general.simplify_full()

    # Apply momentum constraint: D_j Y^j_r = (8πG/c⁴) j_r
    # This gives us the general Einstein tensor component G^t_r

    # From the momentum constraint, we get:
    # DY_general = (8πG/c⁴) e^(-Φ) T^tr
    # Therefore: G^t_r = (c⁴/8πG) DY_general / e^(-Φ) = DY_general * e^Φ * c⁴/(8πG)

    # But we want to express this in terms of ∂_t Φ and ∂_t Λ
    # Substitute the expressions and see what conditions are needed

    return {
        'general_DY': DY_general,
        'general_K_trace': K_trace_general,
        'general_Y_components': (Y_rr_general, Y_thth_general, Y_phph_general),
        'general_K_components': (K_rr_general, K_thth_general, K_phph_general),
        'phi_function': Phi_general,
        'lambda_function': Lambda_general,
        'reduction_condition': 'Λ = -Φ makes e^(2Λ) = e^(-2Φ) = γ_rr^(-1)'
    }

def compute_hamiltonian_constraint(M, A, r, t):
    """
    Compute the Hamiltonian constraint H_⊥ = R^(3) + K² - K_ij K^ij - (16πG/c⁴)ρ = 0
    for our spherical ansatz and verify consistency with flux law
    """
    # Define temporal potential (local variable)
    Phi_local = (1/2) * log(A)

    # Our ansatz: γ_rr = A^(-1), γ_θθ = r², γ_φφ = r²sin²θ
    gamma_rr = A**(-1)
    gamma_thth = r**2
    gamma_phph = r**2 * sin(var('theta'))**2

    # Compute 3D Ricci scalar R^(3) for our metric
    # For diagonal metric γ_ij = diag(A^(-1), r², r²sin²θ)

    # Christoffel symbols for 3-metric
    # Γ^r_rr = (1/2)γ^rr ∂_r γ_rr = (1/2)A(-A^(-2))(-∂_r A) = (1/2A)∂_r A
    Gamma_r_rr = (1/(2*A)) * diff(A, r)

    # Γ^r_θθ = -(1/2)γ^rr ∂_r γ_θθ = -(1/2)A(2r) = -rA
    Gamma_r_thth = -r * A

    # Γ^r_φφ = -(1/2)γ^rr ∂_r γ_φφ = -(1/2)A(2r sin²θ) = -rA sin²θ
    Gamma_r_phph = -r * A * sin(var('theta'))**2

    # Γ^θ_rθ = Γ^θ_θr = (1/2)γ^θθ ∂_r γ_θθ = (1/2)(1/r²)(2r) = 1/r
    Gamma_theta_rtheta = 1/r

    # Γ^φ_rφ = Γ^φ_φr = 1/r
    Gamma_phi_rphi = 1/r

    # Γ^φ_θφ = Γ^φ_φθ = cot(θ)
    Gamma_phi_thetaphi = cos(var('theta'))/sin(var('theta'))

    # 3D Ricci tensor components
    # R_rr = ∂_r Γ^k_rk - ∂_k Γ^k_rr + Γ^k_rk Γ^l_rl - Γ^k_rl Γ^l_rk
    # For spherical symmetry, this simplifies considerably

    # The key insight: for our ansatz A = e^(2Φ), the 3D Ricci scalar is:
    # R^(3) = 2A(∂_r²A/r + ∂_r A/r²) - 2A²(∂_r A)²/r²

    # But we can compute this more systematically using our Φ parameterization
    dPhi_dr = diff(Phi_local, r)
    d2Phi_dr2 = diff(dPhi_dr, r)

    # For metric γ_rr = e^(-2Φ), γ_θθ = r², γ_φφ = r²sin²θ:
    # R^(3) = 4e^(2Φ)[∂_r²Φ + (∂_r Φ)² + (2/r)∂_r Φ]
    R3 = 4 * exp(2*Phi_local) * (d2Phi_dr2 + dPhi_dr**2 + (2/r)*dPhi_dr)

    # Extrinsic curvature components (from our previous calculations)
    dt_Phi = diff(Phi_local, t)
    K_rr = exp(-Phi_local) * exp(-2*Phi_local) * dt_Phi  # = e^(-3Φ) ∂_t Φ
    K_thth = 0  # Time-independent angular components
    K_phph = 0

    # Mixed components
    Krr_up = A * K_rr  # = e^(2Φ) × e^(-3Φ) ∂_t Φ = e^(-Φ) ∂_t Φ
    K_trace = Krr_up  # = e^(-Φ) ∂_t Φ

    # K² term
    K_squared = K_trace**2  # = e^(-2Φ) (∂_t Φ)²

    # K_ij K^ij term
    # K^ij = γ^ik γ^jl K_kl
    # K_ij K^ij = K_rr K^rr + K_θθ K^θθ + K_φφ K^φφ
    # = K_rr × A × K_rr + 0 + 0 = A × K_rr²
    K_ij_Kij = A * K_rr**2
    K_ij_Kij = K_ij_Kij.simplify_full()

    # Hamiltonian constraint (without matter term for now)
    # H_⊥ = R^(3) + K² - K_ij K^ij
    hamiltonian_constraint = R3 + K_squared - K_ij_Kij
    hamiltonian_constraint = hamiltonian_constraint.simplify_full()

    # For vacuum (ρ = 0), this should equal zero
    # For matter, H_⊥ = (16πG/c⁴)ρ

    return {
        'ricci_3d': R3,
        'K_squared': K_squared,
        'K_ij_Kij': K_ij_Kij,
        'hamiltonian_constraint': hamiltonian_constraint,
        'vacuum_condition': 'H_⊥ = 0 gives constraint on A(t,r)',
        'matter_coupling': 'H_⊥ = (16πG/c⁴)ρ for matter',
        'phi_derivatives': (dPhi_dr, d2Phi_dr2, dt_Phi)
    }