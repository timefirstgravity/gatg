#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
ADM (3+1 Decomposition) Equations for GATG
Symbolic representations of ADM decomposition equations
"""

from sage.all import *
from typing import Dict, Any, Optional

class AdmEquations:
    """
    Symbolic equations for ADM 3+1 decomposition in GATG
    """

    def __init__(self, vars_dict: Dict[str, Any]):
        """Initialize with shared variable dictionary"""
        self.vars = vars_dict
        self.equations = {}
        self.category_name = 'adm'

        # Define equations
        self.define_equations()

    def define_equations(self):
        """Define all ADM decomposition equations"""
        r = self.vars['r']
        t = self.vars['t']
        theta = self.vars['theta']

        # General metric forms
        self.equations['metric_general'] = {
            'symbolic': None,  # Line element
            'string': "ds² = -e^(2Φ)dt² + e^(2Λ)dr² + r²dΩ²",
            'latex': r"ds^2 = -e^{2\Phi}dt^2 + e^{2\Lambda}dr^2 + r^2 d\Omega^2",
            'description': "General spherically symmetric metric in exponential form"
        }

        self.equations['metric_ansatz'] = {
            'symbolic': None,  # Line element
            'string': "ds² = -e^(2Φ)dt² + e^(-2Φ)dr² + r²dΩ²",
            'latex': r"ds^2 = -e^{2\Phi}dt^2 + e^{-2\Phi}dr^2 + r^2 d\Omega^2",
            'description': "GATG ansatz metric with Λ = -Φ restriction"
        }

        # Define symbolic functions and tensors
        N = function('N')(t, r)  # Lapse function
        Phi = function('Phi')(t, r)  # Temporal potential
        gamma_rr = function('gamma_rr')(t, r)  # 3-metric component
        gamma_thth = function('gamma_thth')(t, r)
        gamma_phph = function('gamma_phph')(t, r)

        # Define symbolic tensor components
        K_rr = var('K_rr')  # Extrinsic curvature components
        K_thth = var('K_thth')
        K_phph = var('K_phph')

        # Lapse function relation
        self.equations['lapse_function'] = {
            'symbolic': N - exp(Phi),
            'string': "N = e^Φ",
            'latex': r"N = e^\Phi",
            'description': "Lapse function: rate of proper time flow"
        }

        # Normal vector (unit normal to spatial slices)
        self.equations['normal_vector'] = {
            'symbolic': [1/N, 0, 0, 0],  # n^μ components
            'string': "n^μ = (1/N, 0, 0, 0)",
            'latex': r"n^\mu = (1/N, 0, 0, 0)",
            'description': "Unit normal vector to constant-time hypersurfaces"
        }

        # Extrinsic curvature definition
        self.equations['extrinsic_curvature_def'] = {
            'symbolic': K_rr + (1/(2*N)) * diff(gamma_rr, t),
            'string': "K_ij = -(1/2N) ∂_t γ_ij",
            'latex': r"K_{ij} = -\frac{1}{2N} \partial_t \gamma_{ij}",
            'description': "Extrinsic curvature: how spatial slices curve in spacetime"
        }

        # Trace of extrinsic curvature
        K_trace = var('K')
        self.equations['extrinsic_trace'] = {
            'symbolic': K_trace - exp(-Phi) * diff(Phi, t),
            'string': "K = e^(-Φ) ∂_t Φ",
            'latex': r"K = e^{-\Phi} \partial_t \Phi",
            'description': "Trace of extrinsic curvature"
        }

        # 3-metric components (spherical symmetry)
        self.equations['gamma_rr'] = {
            'symbolic': gamma_rr - exp(-2*Phi),
            'string': "γ_rr = e^(-2Φ)",
            'latex': r"\gamma_{rr} = e^{-2\Phi}",
            'description': "Radial component of spatial metric"
        }

        self.equations['gamma_thth'] = {
            'symbolic': r**2,
            'string': "γ_θθ = r²",
            'latex': r"\gamma_{\theta\theta} = r^2",
            'description': "Angular (theta) component of spatial metric"
        }

        self.equations['gamma_phph'] = {
            'symbolic': r**2 * sin(theta)**2,
            'string': "γ_φφ = r²sin²θ",
            'latex': r"\gamma_{\phi\phi} = r^2 \sin^2\theta",
            'description': "Angular (phi) component of spatial metric"
        }

        # Y tensor definition (extrinsic curvature deviation from trace)
        Y_rr = var('Y_rr')
        delta_rr = 1  # Kronecker delta
        self.equations['y_tensor_def'] = {
            'symbolic': Y_rr - (K_rr - delta_rr * K_trace),
            'string': "Y^j_i = K^j_i - δ^j_i K",
            'latex': r"Y^j_i = K^j_i - \delta^j_i K",
            'description': "Y tensor: traceless part of extrinsic curvature"
        }

        # Momentum constraint
        j_r = var('j_r')  # Momentum density
        G = self.vars['G']
        c = self.vars['c']
        self.equations['momentum_constraint'] = {
            'symbolic': var('D_j_Y_jr') - (8*pi*G/c**4) * j_r,
            'string': "D_j Y^j_i = (8πG/c⁴) j_i",
            'latex': r"D_j Y^j_i = \frac{8\pi G}{c^4} j_i",
            'description': "Momentum constraint equation in ADM formulation"
        }

        # Hamiltonian constraint
        R_3 = var('R_3')  # 3D Ricci scalar
        rho = var('rho')  # Energy density
        self.equations['hamiltonian_constraint'] = {
            'symbolic': R_3 + K_trace**2 - var('K_ij_K_ij') - (16*pi*G/c**4)*rho,
            'string': "R^(3) + K² - K_ij K^ij - (16πG/c⁴)ρ = 0",
            'latex': r"R^{(3)} + K^2 - K_{ij} K^{ij} - \frac{16\pi G}{c^4}\rho = 0",
            'description': "Hamiltonian constraint: energy constraint in ADM formulation"
        }

        # ADM line element
        self.equations['adm_line_element'] = {
            'symbolic': None,  # Line element is not a single expression
            'string': "ds² = -N²dt² + γᵢⱼ(dxⁱ + βⁱdt)(dxʲ + βʲdt)",
            'latex': r"ds^2 = -N^2 dt^2 + \gamma_{ij}(dx^i + \beta^i dt)(dx^j + \beta^j dt)",
            'description': "ADM decomposed spacetime line element"
        }

        # Diagonal line element
        self.equations['diagonal_line_element'] = {
            'symbolic': None,  # Line element
            'string': "ds² = -A dt² + A⁻¹ dr² + r² dθ² + r² sin²θ dφ²",
            'latex': r"ds^2 = -A dt^2 + A^{-1} dr^2 + r^2 d\theta^2 + r^2 \sin^2\theta d\phi^2",
            'description': "Diagonal metric line element in spherical coordinates"
        }

        # Gauge choices
        self.equations['zero_shift'] = {
            'symbolic': [0, 0, 0],  # βⁱ = 0
            'string': "βⁱ = 0",
            'latex': r"\beta^i = 0",
            'description': "Zero shift gauge condition"
        }

        self.equations['areal_radius'] = {
            'symbolic': [r**2, r**2 * sin(theta)**2],
            'string': "γ_θθ = r², γ_φφ = r²sin²θ",
            'latex': r"\gamma_{\theta\theta} = r^2, \gamma_{\phi\phi} = r^2\sin^2\theta",
            'description': "Areal radius gauge: sphere area is 4πr²"
        }

        # Ansatz restriction
        Lambda = function('Lambda')(t, r)
        self.equations['ansatz_restriction'] = {
            'symbolic': Lambda + Phi,
            'string': "Λ = -Φ",
            'latex': r"\Lambda = -\Phi",
            'description': "Ansatz restriction relating radial and temporal metric functions"
        }

        # Full ansatz specification
        self.equations['ansatz_full'] = {
            'symbolic': exp(-2*Phi),  # γ_rr = N^(-2) = e^(-2Φ)
            'string': "γ_rr = N^(-2) = e^(-2Φ) [Λ = -Φ restriction]",
            'latex': r"\gamma_{rr} = N^{-2} = e^{-2\Phi} \quad [\Lambda = -\Phi \text{ restriction}]",
            'description': "Complete GATG ansatz specification with restriction noted"
        }

        # Time derivative equations for 3-metric components
        A_function = exp(2*function('Phi')(t, r))
        self.equations['dt_gamma_rr'] = {
            'symbolic': diff(exp(-2*function('Phi')(t, r)), t),
            'string': "∂_t γ_rr = ∂_t(A⁻¹) = ∂_t(e⁻²Φ)",
            'latex': r"\partial_t \gamma_{rr} = \partial_t(A^{-1}) = \partial_t(e^{-2\Phi})",
            'description': "Time derivative of radial metric component"
        }

        self.equations['dt_gamma_theta'] = {
            'symbolic': diff(r**2, t),  # This is 0 since r is spatial
            'string': "∂_t γ_θθ = ∂_t(r²)",
            'latex': r"\partial_t \gamma_{\theta\theta} = \partial_t(r^2)",
            'description': "Time derivative of theta metric component"
        }

        self.equations['dt_gamma_phi'] = {
            'symbolic': diff(r**2 * sin(theta)**2, t),  # This is 0 since r,θ are spatial
            'string': "∂_t γ_φφ = ∂_t(r²sin²θ)",
            'latex': r"\partial_t \gamma_{\phi\phi} = \partial_t(r^2\sin^2\theta)",
            'description': "Time derivative of phi metric component"
        }

        # Extrinsic curvature component equations
        N = exp(function('Phi')(t, r))  # Lapse function
        self.equations['k_rr'] = {
            'symbolic': K_rr + (1/(2*N)) * diff(exp(-2*function('Phi')(t, r)), t),
            'string': "K_rr = -(1/2N) ∂_t γ_rr",
            'latex': r"K_{rr} = -\frac{1}{2N} \partial_t \gamma_{rr}",
            'description': "Extrinsic curvature radial-radial component"
        }

        self.equations['k_theta'] = {
            'symbolic': K_thth + (1/(2*N)) * diff(r**2, t),
            'string': "K_θθ = -(1/2N) ∂_t γ_θθ",
            'latex': r"K_{\theta\theta} = -\frac{1}{2N} \partial_t \gamma_{\theta\theta}",
            'description': "Extrinsic curvature theta-theta component"
        }

        self.equations['k_phi'] = {
            'symbolic': K_phph + (1/(2*N)) * diff(r**2 * sin(theta)**2, t),
            'string': "K_φφ = -(1/2N) ∂_t γ_φφ",
            'latex': r"K_{\phi\phi} = -\frac{1}{2N} \partial_t \gamma_{\phi\phi}",
            'description': "Extrinsic curvature phi-phi component"
        }

        # Y tensor component equations
        Y_rr = var('Y_rr')
        Y_theta = var('Y_theta')
        Y_phi = var('Y_phi')
        self.equations['y_rr_component'] = {
            'symbolic': Y_rr - (K_rr - K_trace),
            'string': "Y^r_r = K^r_r - δ^r_r K = K^r_r - K",
            'latex': r"Y^r_r = K^r_r - \delta^r_r K = K^r_r - K",
            'description': "Y tensor radial-radial component"
        }

        self.equations['y_theta_component'] = {
            'symbolic': Y_theta + K_trace,  # Y^θ_θ = 0 - K = -K
            'string': "Y^θ_θ = K^θ_θ - δ^θ_θ K = 0 - K",
            'latex': r"Y^\theta_\theta = K^\theta_\theta - \delta^\theta_\theta K = 0 - K",
            'description': "Y tensor theta-theta component"
        }

        self.equations['y_phi_component'] = {
            'symbolic': Y_phi + K_trace,  # Y^φ_φ = 0 - K = -K
            'string': "Y^φ_φ = K^φ_φ - δ^φ_φ K = 0 - K",
            'latex': r"Y^\phi_\phi = K^\phi_\phi - \delta^\phi_\phi K = 0 - K",
            'description': "Y tensor phi-phi component"
        }

        # Additional Y tensor relations
        self.equations['y_traceless_relation'] = {
            'symbolic': None,  # General traceless tensor definition
            'string': "Ỹ^j_i = K^j_i - (1/3)δ^j_i K",
            'latex': r"\tilde{Y}^j_i = K^j_i - \frac{1}{3}\delta^j_i K",
            'description': "Alternative traceless decomposition of extrinsic curvature"
        }

        # Temporal potential definition
        self.equations['temporal_potential'] = {
            'symbolic': None,  # Definition
            'string': "A ≡ e^(2Φ)",
            'latex': r"A \equiv e^{2\Phi}",
            'description': "Temporal potential A in terms of lapse function exponent"
        }

        # A function definition (alternative form)
        self.equations['a_function_definition'] = {
            'symbolic': None,  # Definition
            'string': "A = e^(2Φ)",
            'latex': r"A = e^{2\Phi}",
            'description': "A function defined as exponential of twice the temporal potential"
        }

        # Hamiltonian constraints
        rho = var('rho', domain='positive')
        self.equations['vacuum_hamiltonian'] = {
            'symbolic': None,  # Constraint equation
            'string': "H_⊥ = 0",
            'latex': r"H_\perp = 0",
            'description': "Vacuum Hamiltonian constraint (no matter)"
        }

        self.equations['matter_hamiltonian'] = {
            'symbolic': (16*pi*G/c**4) * rho,
            'string': "H_⊥ = (16πG/c⁴)ρ",
            'latex': r"H_\perp = \frac{16\pi G}{c^4}\rho",
            'description': "Matter Hamiltonian constraint with energy density"
        }

        # Mathematical definitions
        self.equations['kronecker_delta'] = {
            'symbolic': None,  # Mathematical definition
            'string': "δ^j_i (= 1 if i=j, = 0 if i≠j)",
            'latex': r"\delta^j_i \text{ (= 1 if } i=j\text{, = 0 if } i \neq j\text{)}",
            'description': "Kronecker delta tensor definition"
        }

        # Index raising formula
        self.equations['index_raising_formula'] = {
            'symbolic': None,  # General formula
            'string': "K^i_j = γ^ik K_kj",
            'latex': r"K^i_j = \gamma^{ik} K_{kj}",
            'description': "Index raising formula using inverse 3-metric"
        }

        # Inverse 3-metric components
        self.equations['gamma_inv_rr'] = {
            'symbolic': exp(2*Phi),  # A = e^(2Φ)
            'string': "γ^rr = A = e^(2Φ)",
            'latex': r"\gamma^{rr} = A = e^{2\Phi}",
            'description': "Inverse radial component of 3-metric"
        }

        self.equations['gamma_inv_theta'] = {
            'symbolic': 1/r**2,
            'string': "γ^θθ = 1/r²",
            'latex': r"\gamma^{\theta\theta} = \frac{1}{r^2}",
            'description': "Inverse theta component of 3-metric"
        }

        self.equations['gamma_inv_phi'] = {
            'symbolic': 1/(r**2 * sin(theta)**2),
            'string': "γ^φφ = 1/(r²sin²θ)",
            'latex': r"\gamma^{\phi\phi} = \frac{1}{r^2\sin^2\theta}",
            'description': "Inverse phi component of 3-metric"
        }

        # Mixed extrinsic curvature components
        self.equations['k_mixed_rr'] = {
            'symbolic': None,  # γ^rr K_rr
            'string': "K^r_r = γ^rr K_rr",
            'latex': r"K^r_r = \gamma^{rr} K_{rr}",
            'description': "Mixed radial-radial extrinsic curvature component"
        }

        self.equations['k_mixed_theta'] = {
            'symbolic': None,  # γ^θθ K_θθ
            'string': "K^θ_θ = γ^θθ K_θθ",
            'latex': r"K^\theta_\theta = \gamma^{\theta\theta} K_{\theta\theta}",
            'description': "Mixed theta-theta extrinsic curvature component"
        }

        self.equations['k_mixed_phi'] = {
            'symbolic': None,  # γ^φφ K_φφ
            'string': "K^φ_φ = γ^φφ K_φφ",
            'latex': r"K^\phi_\phi = \gamma^{\phi\phi} K_{\phi\phi}",
            'description': "Mixed phi-phi extrinsic curvature component"
        }

        # Extrinsic curvature trace definition
        self.equations['extrinsic_trace_definition'] = {
            'symbolic': None,  # Sum of mixed components
            'string': "K ≡ K^i_i = K^r_r + K^θ_θ + K^φ_φ",
            'latex': r"K \equiv K^i_i = K^r_r + K^\theta_\theta + K^\phi_\phi",
            'description': "Trace of extrinsic curvature as sum of mixed diagonal components"
        }

        # Trace formula verification
        self.equations['trace_formula_verification'] = {
            'symbolic': K_trace == exp(-Phi) * diff(Phi, t),
            'string': "K = e^(-Φ) ∂_t Φ",
            'latex': r"K = e^{-\Phi} \partial_t \Phi",
            'description': "Trace formula verification: relation between trace and temporal potential evolution"
        }

    def get_equations(self) -> Dict[str, Dict[str, Any]]:
        """Return all equations in this category"""
        return self.equations.copy()

    def get_equation_keys(self) -> list:
        """Return list of equation keys in this category"""
        return list(self.equations.keys())

# Export the class
__all__ = ['AdmEquations']