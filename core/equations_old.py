#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Standard Equations
Centralized definitions of all mathematical equations used throughout GATG modules
This eliminates duplication and ensures consistency across all print functions

NOTE: This module now integrates with symbolic_equations.py for equations that have
been migrated to the symbolic framework. This provides backward compatibility while
enabling mathematical verification of documented equations.
"""

# Try to import symbolic equations for enhanced functionality
try:
    # Try relative import first (when used as module)
    from .symbolic_equations import symbolic_equations, SymbolicGATGEquations
    SYMBOLIC_AVAILABLE = True
except ImportError:
    try:
        # Try absolute import (when run directly)
        from symbolic_equations import symbolic_equations, SymbolicGATGEquations
        SYMBOLIC_AVAILABLE = True
    except ImportError:
        SYMBOLIC_AVAILABLE = False

class GATGEquations:
    """Centralized repository of all GATG equations as formatted strings"""

    # Note: Most descriptive text equations have been inlined directly into print functions
    # where they are used, as they are not mathematical equations to be computed.
    # This centralized equations module now focuses on equations that are:
    # 1. Used in multiple places across modules
    # 2. Mathematical equations that can be symbolically verified
    # 3. Complex enough to warrant centralization

    @classmethod
    def get_flux_law_boxed(cls):
        """Return the flux law in a formatted box"""
        return f"""┌{"─"*50}┐
│{cls.FLUX_LAW:^50}│
└{"─"*50}┘"""

    @classmethod
    def get_metric_forms(cls):
        """Return all metric forms as a formatted list"""
        return f"""General: {cls.METRIC_GENERAL}
Ansatz:  {cls.METRIC_ANSATZ}
Vaidya (EF): {cls.METRIC_VAIDYA_EF}
Vaidya (diagonal): {cls.METRIC_VAIDYA_DIAGONAL}"""

    # get_adm_variables() method removed - ADM equations now available in symbolic_equations.categories['adm']

    @classmethod
    def get_symbolic(cls, equation_name):
        """
        Get symbolic representation of an equation if available

        This provides backward compatibility with the new symbolic framework.
        Returns the symbolic expression if the equation has been migrated,
        otherwise returns None.
        """
        if not SYMBOLIC_AVAILABLE:
            return None

        # Map class variable names to symbolic equation keys
        equation_map = {
            # Schwarzschild equations
            'SCHWARZSCHILD_ODE': 'schwarzschild_ode',
            'SCHWARZSCHILD_CONSTRAINT': 'schwarzschild_constraint',
            'SCHWARZSCHILD_SOLUTION': 'schwarzschild_solution',
            'SCHWARZSCHILD_RADIUS': 'schwarzschild_radius',
            'SCHWARZSCHILD_METRIC': 'schwarzschild_metric',
            'SCHWARZSCHILD_METRIC_GENERAL': 'schwarzschild_metric_general',
            'SCHWARZSCHILD_RECOVERY': 'schwarzschild_recovery',
            'SCHWARZSCHILD_FORM': 'schwarzschild_form',
            'VACUUM_LIMIT': 'vacuum_limit',
            # Flux law equations
            'FLUX_LAW': 'flux_law',
            'FLUX_LAW_GEOMETRIZED': 'flux_law_geometrized',
            'EINSTEIN_TR': 'einstein_tr',
            'EINSTEIN_CONSTRAINT': 'einstein_constraint',
            # General metric forms
            'METRIC_GENERAL': 'metric_general',
            'METRIC_ANSATZ': 'metric_ansatz',
            'METRIC_VAIDYA_EF': 'vaidya_metric_ingoing',
            'METRIC_VAIDYA_DIAGONAL': 'vaidya_metric_diagonal',
            # ADM equations
            'LAPSE_FUNCTION': 'lapse_function',
            'NORMAL_VECTOR': 'normal_vector',
            'EXTRINSIC_CURVATURE_DEF': 'extrinsic_curvature_def',
            'EXTRINSIC_TRACE': 'extrinsic_trace',
            'GAMMA_RR': 'gamma_rr',
            'GAMMA_THTH': 'gamma_thth',
            'GAMMA_PHPH': 'gamma_phph',
            'Y_TENSOR_DEF': 'y_tensor_def',
            'Y_RR_ZERO': 'y_rr_component',
            'Y_TRACELESS_RELATION': 'y_traceless_relation',
            'TEMPORAL_POTENTIAL': 'temporal_potential',
            'MOMENTUM_CONSTRAINT': 'momentum_constraint',
            'HAMILTONIAN_CONSTRAINT': 'hamiltonian_constraint',
            'VACUUM_HAMILTONIAN': 'vacuum_hamiltonian',
            'MATTER_HAMILTONIAN': 'matter_hamiltonian',
            'KRONECKER_DELTA': 'kronecker_delta',
            'ZERO_SHIFT': 'zero_shift',
            'AREAL_RADIUS': 'areal_radius',
            'ANSATZ_RESTRICTION': 'ansatz_restriction',
            'ADM_GENERAL_LINE_ELEMENT': 'adm_line_element',
            'ADM_LINE_ELEMENT_DETAILED': 'adm_line_element',
            'DIAGONAL_LINE_ELEMENT': 'diagonal_line_element',
            'MIXED_COMPONENT_FORMULA': 'index_raising_formula',
            'TRACE_FORMULA_VERIFICATION': 'trace_formula_verification',
            'DT_GAMMA_RR_EQUATION': 'dt_gamma_rr',
            'DT_GAMMA_THETA_EQUATION': 'dt_gamma_theta',
            'DT_GAMMA_PHI_EQUATION': 'dt_gamma_phi',
            'K_RR_EQUATION': 'k_rr',
            'K_THETA_EQUATION': 'k_theta',
            'K_PHI_EQUATION': 'k_phi',
            'INDEX_RAISING_FORMULA': 'index_raising_formula',
            'GAMMA_INV_RR': 'gamma_inv_rr',
            'GAMMA_INV_THETA': 'gamma_inv_theta',
            'GAMMA_INV_PHI': 'gamma_inv_phi',
            'K_MIXED_RR': 'k_mixed_rr',
            'K_MIXED_THETA': 'k_mixed_theta',
            'K_MIXED_PHI': 'k_mixed_phi',
            'EXTRINSIC_TRACE_DEFINITION': 'extrinsic_trace_definition',
            'Y_RR_COMPONENT': 'y_rr_component',
            'Y_THETA_COMPONENT': 'y_theta_component',
            'Y_PHI_COMPONENT': 'y_phi_component',
            # Christoffel symbols
            'CHRISTOFFEL_T_TT': 'christoffel_t_tt',
            'CHRISTOFFEL_T_TR': 'christoffel_t_tr',
            'CHRISTOFFEL_T_RR': 'christoffel_t_rr',
            'CHRISTOFFEL_R_TT': 'christoffel_r_tt',
            'CHRISTOFFEL_R_TR': 'christoffel_r_tr',
            'CHRISTOFFEL_R_RR': 'christoffel_r_rr',
            'CHRISTOFFEL_THETA_RTHETA': 'christoffel_theta_rtheta',
            'CHRISTOFFEL_PHI_RPHI': 'christoffel_phi_rphi',
            'CHRISTOFFEL_R_THETATHETA': 'christoffel_r_thetatheta',
            'CHRISTOFFEL_R_PHIPHI': 'christoffel_r_phiphi',
            'RICCI_3D_SPHERICAL': 'ricci_3d_spherical',
            # Kerr equations
            'KERR_SIGMA': 'kerr_sigma',
            'KERR_DELTA': 'kerr_delta',
            'KERR_G_TT': 'kerr_g_tt',
            'KERR_G_RR': 'kerr_g_rr',
            'KERR_G_THTH': 'kerr_g_thth',
            'KERR_G_PHPH': 'kerr_g_phph',
            'KERR_G_TPH': 'kerr_g_tph',
            'KERR_EVENT_HORIZONS': 'kerr_event_horizons',
            'KERR_ERGOSPHERE': 'kerr_ergosphere',
            'KERR_ANGULAR_MOMENTUM': 'kerr_angular_momentum',
            'KERR_SCHWARZSCHILD_LIMIT': 'kerr_schwarzschild_limit',
            'KERR_XI': 'kerr_xi',
            'KERR_LAPSE_FUNCTION': 'kerr_lapse_function',
            'KERR_LAPSE_SIMPLIFIED': 'kerr_lapse_simplified',
            'KERR_SHIFT_VECTOR': 'kerr_shift_vector',
            # Cosmological equations
            'FRW_METRIC': 'frw_metric',
            'FRW_METRIC_FLAT': 'frw_metric_flat',
            'FRW_LAPSE_FIRST': 'frw_lapse_first',
            'HUBBLE_PARAMETER': 'hubble_parameter',
            'SCALE_FACTOR_LAPSE_FIRST': 'scale_factor_lapse_first',
            'HUBBLE_LAPSE_FIRST': 'hubble_lapse_first',
            'FRIEDMANN_FIRST': 'friedmann_first',
            'FRIEDMANN_SECOND': 'friedmann_second',
            'FRIEDMANN_FIRST_NORMALIZED': 'friedmann_first_normalized',
            'MATTER_DOMINATED': 'matter_dominated',
            'RADIATION_DOMINATED': 'radiation_dominated',
            'LAMBDA_CDM': 'lambda_cdm',
            'REDSHIFT_SCALE_RELATION': 'redshift_scale_relation',
            'EZ_LAMBDA_CDM': 'ez_lambda_cdm',
            'LUMINOSITY_DISTANCE': 'luminosity_distance',
            'ANGULAR_DIAMETER_DISTANCE': 'angular_diameter_distance',
            'AGE_OF_UNIVERSE': 'age_of_universe',
            # Conservation equations
            'CONSERVATION_EXPANDED': 'conservation_expanded',
            'CONSERVATION_EXPLICIT': 'conservation_explicit',
            'STRESS_ENERGY_CONSERVATION': 'stress_energy_conservation',
            'COVARIANT_DIVERGENCE': 'covariant_divergence',
            'MIXED_STRESS_ENERGY': 'mixed_stress_energy',
            'MIXED_EINSTEIN': 'mixed_einstein',
            'COVARIANT_DIVERGENCE_FORMULA': 'covariant_divergence_formula',
            'COVARIANT_DIVERGENCE_SIMPLIFIED': 'covariant_divergence_simplified',
            'MOMENTUM_CONSTRAINT_PRECISE': 'momentum_constraint_precise',
            'MOMENTUM_CONSTRAINT_DETAILED': 'momentum_constraint_detailed',
            'EXTRINSIC_CANCELLATION': 'extrinsic_cancellation',
            'GENERAL_EINSTEIN_TR': 'general_einstein_tr',
            'GENERAL_EINSTEIN_TR_FORMULA': 'general_einstein_tr',
            'GENERAL_EINSTEIN_REDUCTION': 'general_einstein_reduction',
            'DIVERGENCE_SIMPLIFIED_FORMULA': 'divergence_simplified_formula',
            'DIVERGENCE_APPLIED_FORMULA': 'divergence_applied_formula',
            'MOMENTUM_CONSTRAINT_MATHEMATICAL': 'momentum_constraint_precise',
            'J_R_PROJECTION_STEP1': 'j_r_projection_step1',
            'J_R_PROJECTION_STEP2': 'j_r_projection_step2',
            'J_R_PROJECTION_RESULT': 'j_r_projection_result',
            'DIVERGENCE_RESULT': 'divergence_result',
            'CONSTRAINT_REQUIREMENT': 'constraint_requirement',
            'CONSTRAINT_EQUALITY': 'constraint_equality',
            'MOMENTUM_CONSTRAINT_DETAILED_LONG': 'momentum_constraint_detailed',
            'CANCELLATION_RESULT': 'cancellation_result',
            'EINSTEIN_CONSTRAINT_EQUATION': 'einstein_constraint',
            # Coordinate transformation equations
            'VAIDYA_TRANSFORMATION': 'vaidya_transformation',
            'VAIDYA_INTEGRATION': 'vaidya_integration',
            'COORDINATE_JACOBIAN': 'coordinate_jacobian',
            'EF_TO_DIAGONAL_CONDITION': 'ef_to_diagonal_condition',
            'COORDINATE_REDEFINITION': 'coordinate_redefinition',
            'TIME_REDEFINITION': 'time_redefinition',
            'VAIDYA_A_FUNCTION': 'vaidya_a_function',
            'VAIDYA_PHI_RELATION': 'vaidya_phi_relation',
            'PAINLEVE_GULLSTRAND_METRIC': 'painleve_gullstrand_metric',
            'PAINLEVE_GULLSTRAND_VELOCITY': 'painleve_gullstrand_velocity',
            'EDDINGTON_FINKELSTEIN_INGOING': 'eddington_finkelstein_ingoing',
            'EDDINGTON_FINKELSTEIN_OUTGOING': 'eddington_finkelstein_outgoing',
            'TORTOISE_COORDINATE': 'tortoise_coordinate',
            'KRUSKAL_SZEKERES_U': 'kruskal_szekeres_u',
            'KRUSKAL_SZEKERES_V': 'kruskal_szekeres_v',
            'KRUSKAL_SZEKERES_METRIC': 'kruskal_szekeres_metric',
            # Reissner-Nordström equations
            'REISSNER_NORDSTROM_METRIC': 'reissner_nordstrom_metric',
            'REISSNER_NORDSTROM_SOLUTION': 'reissner_nordstrom_solution',
            'REISSNER_NORDSTROM_HORIZONS': 'reissner_nordstrom_horizons',
            'REISSNER_NORDSTROM_CHARGE': 'reissner_nordstrom_charge',
            # Vaidya equations
            'VAIDYA_METRIC_INGOING': 'vaidya_metric_ingoing',
            'VAIDYA_METRIC_OUTGOING': 'vaidya_metric_outgoing',
            'VAIDYA_MASS_FUNCTION': 'vaidya_mass_function',
            'VAIDYA_A_FUNCTION': 'vaidya_a_function',
            'VAIDYA_PHI_RELATION': 'vaidya_phi_relation',
            'MASS_FUNCTION_LINEAR': 'mass_function_linear',
            'A_FUNCTION_DEFINITION': 'a_function_definition',
            # Kerr-Newman equations
            'KERR_NEWMAN_DELTA': 'kerr_newman_delta',
            'KERR_NEWMAN_HORIZONS': 'kerr_newman_horizons',
            # CMB observable equations
            'SOUND_HORIZON': 'sound_horizon',
            'THETA_STAR': 'theta_star',
            'SAHA_EQUATION': 'saha_equation',
            'AGE_OF_UNIVERSE': 'age_of_universe',
            'ANGULAR_DIAMETER_DISTANCE': 'angular_diameter_distance',
        }

        if equation_name in equation_map:
            try:
                return symbolic_equations.get_symbolic(equation_map[equation_name])
            except KeyError:
                return None

        return None

    @classmethod
    def verify_equation(cls, equation_name):
        """
        Verify an equation using the symbolic framework if available

        Returns:
            bool: True if equation is verified, False if not verified,
                  None if symbolic framework not available
        """
        if not SYMBOLIC_AVAILABLE:
            return None

        # Specific verification methods
        if equation_name == 'SCHWARZSCHILD_ODE':
            return symbolic_equations.verify_schwarzschild_ode()
        elif equation_name == 'FLUX_LAW':
            return symbolic_equations.verify_flux_law_consistency()
        elif equation_name == 'KERR_SCHWARZSCHILD_LIMIT':
            return symbolic_equations.verify_kerr_schwarzschild_limit()
        elif equation_name == 'KERR_EVENT_HORIZONS':
            return symbolic_equations.verify_kerr_horizons()
        elif equation_name == 'KERR_SIGMA' or equation_name == 'KERR_DELTA':
            return symbolic_equations.verify_kerr_metric_consistency()
        elif equation_name == 'HUBBLE_LAPSE_FIRST' or equation_name == 'SCALE_FACTOR_LAPSE_FIRST':
            return symbolic_equations.verify_lapse_first_hubble_relation()
        elif equation_name == 'MATTER_DOMINATED':
            return symbolic_equations.verify_matter_scaling()
        elif equation_name == 'RADIATION_DOMINATED':
            return symbolic_equations.verify_radiation_scaling()
        elif equation_name == 'REDSHIFT_SCALE_RELATION':
            return symbolic_equations.verify_cosmological_redshift_relation()
        # Christoffel symbol verifications
        elif equation_name.startswith('CHRISTOFFEL_'):
            # All Christoffel symbols can be verified using symmetry checks
            return symbolic_equations.verify_christoffel_spherical_symmetry()

        return None

    @classmethod
    def has_symbolic(cls, equation_name):
        """Check if an equation has a symbolic representation available"""
        return cls.get_symbolic(equation_name) is not None

# Create module-level constants for easy import
FLUX_LAW = GATGEquations.FLUX_LAW
EINSTEIN_TR = GATGEquations.EINSTEIN_TR
METRIC_ANSATZ = GATGEquations.METRIC_ANSATZ
EXTRINSIC_CURVATURE_DEF = GATGEquations.EXTRINSIC_CURVATURE_DEF
MOMENTUM_CONSTRAINT = GATGEquations.MOMENTUM_CONSTRAINT

# Cosmological constants
FRW_METRIC = GATGEquations.FRW_METRIC
FRIEDMANN_FIRST = GATGEquations.FRIEDMANN_FIRST
FRIEDMANN_SECOND = GATGEquations.FRIEDMANN_SECOND
SCALE_FACTOR_LAPSE_FIRST = GATGEquations.SCALE_FACTOR_LAPSE_FIRST
HUBBLE_LAPSE_FIRST = GATGEquations.HUBBLE_LAPSE_FIRST

# Additional module-level constants for new metrics
SCHWARZSCHILD_METRIC = GATGEquations.SCHWARZSCHILD_METRIC
SCHWARZSCHILD_ODE = GATGEquations.SCHWARZSCHILD_ODE
KERR_METRIC = GATGEquations.KERR_METRIC
KERR_SIGMA = GATGEquations.KERR_SIGMA
KERR_DELTA = GATGEquations.KERR_DELTA
VAIDYA_METRIC_INGOING = GATGEquations.VAIDYA_METRIC_INGOING
REISSNER_NORDSTROM_METRIC = GATGEquations.REISSNER_NORDSTROM_METRIC

# Export the class and key constants
__all__ = [
    'GATGEquations',
    'FLUX_LAW',
    'EINSTEIN_TR',
    'METRIC_ANSATZ',
    'EXTRINSIC_CURVATURE_DEF',
    'MOMENTUM_CONSTRAINT',
    'FRW_METRIC',
    'FRIEDMANN_FIRST',
    'FRIEDMANN_SECOND',
    'SCALE_FACTOR_LAPSE_FIRST',
    'HUBBLE_LAPSE_FIRST',
    'SCHWARZSCHILD_METRIC',
    'SCHWARZSCHILD_ODE',
    'KERR_METRIC',
    'KERR_SIGMA',
    'KERR_DELTA',
    'VAIDYA_METRIC_INGOING',
    'REISSNER_NORDSTROM_METRIC'
]
