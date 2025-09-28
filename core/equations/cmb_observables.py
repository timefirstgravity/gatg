#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
CMB and Cosmological Observables Equations for GATG
Symbolic representations of cosmic microwave background and observational cosmology
"""

from sage.all import *
from typing import Dict, Any, Optional

class CmbObservablesEquations:
    """
    Symbolic equations for CMB physics and cosmological observables in GATG
    """

    def __init__(self, vars_dict: Dict[str, Any]):
        """Initialize with shared variable dictionary"""
        self.vars = vars_dict
        self.equations = {}
        self.category_name = 'cmb_observables'

        # Define equations
        self.define_equations()

    def define_equations(self):
        """Define CMB and observational cosmology equations"""

        # Additional variables for CMB physics
        z_star = var('z_star', domain='positive')  # Recombination redshift
        z = self.vars['z']  # Redshift (already defined)
        H_0 = self.vars['H_0']  # Hubble constant
        c = self.vars['c']  # Speed of light

        # CMB-specific variables
        c_s = var('c_s', domain='positive')  # Sound speed
        eta = var('eta', domain='real')  # Conformal time
        r_s = var('r_s', domain='positive')  # Sound horizon
        theta_star = var('theta_star', domain='positive')  # Angular scale at recombination

        # Recombination physics variables
        x_e = var('x_e', domain='positive')  # Ionization fraction
        n_H = var('n_H', domain='positive')  # Hydrogen density
        T = var('T', domain='positive')  # Temperature
        m_e = var('m_e', domain='positive')  # Electron mass
        k = var('k', domain='positive')  # Boltzmann constant
        h = var('h', domain='positive')  # Planck constant

        # Sound horizon at recombination
        self.equations['sound_horizon'] = {
            'symbolic': None,  # Integral expression
            'string': "r_s = ∫₀^z_* c_s dη",
            'latex': r"r_s = \int_0^{z_*} c_s \, d\eta",
            'description': "Sound horizon: distance sound waves could travel by recombination"
        }

        # Angular size of sound horizon
        self.equations['theta_star'] = {
            'symbolic': r_s / var('d_A'),  # r_s/d_A(z_*)
            'string': "θ_* = r_s(z_*)/d_A(z_*)",
            'latex': r"\theta_* = \frac{r_s(z_*)}{d_A(z_*)}",
            'description': "Angular scale of sound horizon at recombination (first acoustic peak)"
        }

        # Saha equation for recombination
        self.equations['saha_equation'] = {
            'symbolic': x_e**2/(1-x_e) - (2*pi*m_e*k*T/h**2)**(3/2) * exp(-13.6*var('eV')/(k*T)) / n_H,
            'string': "x_e²/(1-x_e) = (2πm_e kT/h²)^(3/2) e^(-13.6eV/kT) / n_H",
            'latex': r"\frac{x_e^2}{1-x_e} = \left(\frac{2\pi m_e kT}{h^2}\right)^{3/2} \frac{e^{-13.6\text{eV}/kT}}{n_H}",
            'description': "Saha equation governing hydrogen recombination"
        }

        # Recombination redshift (approximate)
        self.equations['recombination_redshift'] = {
            'symbolic': 1100,  # Approximate value
            'string': "z_* ≈ 1100",
            'latex': r"z_* \approx 1100",
            'description': "Redshift of recombination (when CMB was last scattered)"
        }

        # CMB temperature today
        self.equations['cmb_temperature_today'] = {
            'symbolic': 2.725,  # Kelvin
            'string': "T_CMB = 2.725 K",
            'latex': r"T_{\text{CMB}} = 2.725 \text{ K}",
            'description': "Current CMB temperature (COBE/WMAP/Planck measurement)"
        }

        # CMB temperature scaling
        self.equations['cmb_temperature_scaling'] = {
            'symbolic': var('T_CMB') * (1 + z),
            'string': "T(z) = T_CMB × (1 + z)",
            'latex': r"T(z) = T_{\text{CMB}} \times (1 + z)",
            'description': "CMB temperature as function of redshift"
        }

        # Acoustic oscillations
        self.equations['acoustic_oscillation_scale'] = {
            'symbolic': pi / theta_star,
            'string': "l_acoustic ≈ π/θ_*",
            'latex': r"l_{\text{acoustic}} \approx \frac{\pi}{\theta_*}",
            'description': "Multipole scale of first acoustic peak in CMB"
        }

        # Silk damping scale
        self.equations['silk_damping_scale'] = {
            'symbolic': None,  # Complex expression
            'string': "l_Silk ∝ (η_* × k_Silk)^(-1)",
            'latex': r"l_{\text{Silk}} \propto (\eta_* \times k_{\text{Silk}})^{-1}",
            'description': "Silk damping scale (small-scale suppression in CMB)"
        }

        # Distance measures in cosmology
        self.equations['comoving_distance'] = {
            'symbolic': None,  # Integral
            'string': "d_C(z) = ∫₀^z c dz'/H(z')",
            'latex': r"d_C(z) = \int_0^z \frac{c \, dz'}{H(z')}",
            'description': "Comoving distance to redshift z"
        }

        self.equations['angular_diameter_distance'] = {
            'symbolic': None,  # From cosmological module, but repeated for completeness
            'string': "d_A(z) = d_L(z)/(1+z)²",
            'latex': r"d_A(z) = \frac{d_L(z)}{(1+z)^2}",
            'description': "Angular diameter distance (for CMB angular scales)"
        }

        # Baryon acoustic oscillations (BAO)
        self.equations['bao_standard_ruler'] = {
            'symbolic': r_s,  # Same as sound horizon
            'string': "r_BAO = r_s(z_drag)",
            'latex': r"r_{\text{BAO}} = r_s(z_{\text{drag}})",
            'description': "BAO standard ruler (sound horizon at drag epoch)"
        }

        # CMB lensing
        self.equations['lensing_potential'] = {
            'symbolic': None,  # Complex integral
            'string': "φ(n̂) = -2∫₀^η_* dη (η_*-η)/η_* × Ψ(η×n̂,η)",
            'latex': r"\varphi(\hat{n}) = -2\int_0^{\eta_*} d\eta \frac{\eta_* - \eta}{\eta_*} \times \Psi(\eta \hat{n}, \eta)",
            'description': "CMB lensing potential from matter fluctuations"
        }

        # Power spectrum definitions
        self.equations['temperature_power_spectrum'] = {
            'symbolic': None,  # Statistical quantity
            'string': "C_l^TT = ⟨a_lm^T a_l'm'*⟩ δ_ll' δ_mm'",
            'latex': r"C_l^{TT} = \langle a_{lm}^T a_{l'm'}^* \rangle \delta_{ll'} \delta_{mm'}",
            'description': "CMB temperature angular power spectrum"
        }

        self.equations['polarization_power_spectrum'] = {
            'symbolic': None,  # Statistical quantity
            'string': "C_l^EE = ⟨a_lm^E a_l'm'*⟩ δ_ll' δ_mm'",
            'latex': r"C_l^{EE} = \langle a_{lm}^E a_{l'm'}^* \rangle \delta_{ll'} \delta_{mm'}",
            'description': "CMB E-mode polarization power spectrum"
        }

        # Cosmological parameter constraints
        self.equations['cmb_parameter_constraints'] = {
            'symbolic': None,  # Observational results
            'string': "Ω_b h² = 0.02237, Ω_c h² = 0.1200, H₀ = 67.4 km/s/Mpc",
            'latex': r"\Omega_b h^2 = 0.02237, \quad \Omega_c h^2 = 0.1200, \quad H_0 = 67.4 \text{ km/s/Mpc}",
            'description': "Key cosmological parameters from Planck CMB observations"
        }

        # Alternative dark energy models (temporal dark energy from GATG)
        self.equations['temporal_dark_energy_ez'] = {
            'symbolic': None,  # Model-dependent
            'string': "E²(z) = Ω_m(1+z)³ + ρ_TDE",
            'latex': r"E^2(z) = \Omega_m(1+z)^3 + \rho_{\text{TDE}}",
            'description': "E(z) function for temporal dark energy models"
        }

        self.equations['temporal_dark_energy_density'] = {
            'symbolic': None,  # Alternative to Λ
            'string': "ρ = ρ_m a^(-3) + ρ_TDE",
            'latex': r"\rho = \rho_m a^{-3} + \rho_{\text{TDE}}",
            'description': "Energy density with temporal dark energy instead of cosmological constant"
        }

        # Observational tests
        self.equations['supernovae_distance_modulus'] = {
            'symbolic': 5 * log(var('d_L'), 10) + 25,  # Distance modulus
            'string': "μ(z) = 5 log₁₀(d_L/Mpc) + 25",
            'latex': r"\mu(z) = 5 \log_{10}(d_L/\text{Mpc}) + 25",
            'description': "Distance modulus for Type Ia supernovae"
        }

        self.equations['age_constraint'] = {
            'symbolic': None,  # Integral constraint
            'string': "t₀ = ∫₀^∞ dz/[(1+z)H₀E(z)] > 13.8 Gyr",
            'latex': r"t_0 = \int_0^\infty \frac{dz}{(1+z)H_0 E(z)} > 13.8 \text{ Gyr}",
            'description': "Age of universe constraint from stellar evolution"
        }

    def verify_recombination_physics(self) -> bool:
        """Verify basic recombination physics relationships"""
        # Recombination should occur around z ~ 1100
        # CMB temperature should scale as (1+z)

        recomb_z = self.equations['recombination_redshift']['symbolic']
        temp_scaling = self.equations['cmb_temperature_scaling']['string']

        # Basic sanity checks
        return recomb_z > 1000 and "1 + z" in temp_scaling

    def verify_acoustic_scale_consistency(self) -> bool:
        """Verify that acoustic scales are consistent"""
        # The acoustic scale should be related to sound horizon and angular diameter distance

        theta_def = self.equations['theta_star']['string']
        acoustic_scale = self.equations['acoustic_oscillation_scale']['string']

        # Should involve sound horizon and angular diameter distance
        return "r_s" in theta_def and "d_A" in theta_def and "θ_*" in acoustic_scale

    def verify_distance_measure_relationships(self) -> bool:
        """Verify relationships between cosmological distance measures"""
        # Angular diameter distance should be related to luminosity distance

        d_A_relation = self.equations['angular_diameter_distance']['string']

        # Should have the (1+z)² relationship
        return "d_L(z)/(1+z)²" in d_A_relation

    def get_equations(self) -> Dict[str, Dict[str, Any]]:
        """Return all equations in this category"""
        return self.equations.copy()

    def get_equation_keys(self) -> list:
        """Return list of equation keys in this category"""
        return list(self.equations.keys())

# Export the class
__all__ = ['CmbObservablesEquations']