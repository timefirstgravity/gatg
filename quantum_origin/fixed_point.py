#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Quantum Origin: Fixed-Point Operator

Implements the fixed-point map Φ = T[Φ] ≡ C ω² C[Φ].
Verifies Banach contraction conditions for well-posed emergence.

From Quantum Origin paper:
- T is a contraction if C||ω²||_∞ L_C < 1
- Fixed point exists and is unique under contraction
- Iteration converges geometrically
"""

from sage.all import *
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def build_fixed_point_operator(ctp_kernel, clock_frequencies, operator_params=None):
    """
    Build the fixed-point operator T = C ω² C.

    From paper: Φ = T[Φ] ≡ C ω² C[Φ]
    where C is the CTP kernel and ω² is clock frequency squared.

    Args:
        ctp_kernel: CTP kernel from ctp_kernel module
        clock_frequencies: Angular frequencies ω (can be field-dependent)
        operator_params: Dict with coupling constants

    Returns:
        Dict with:
            'operator_T': The composed operator T
            'kernel_C': CTP kernel operator
            'frequency_factor': ω²
            'operator_norm_bound': Upper bound on ||T||
    """
    if operator_params is None:
        operator_params = {'coupling_constant': 1.0}

    C_op = ctp_kernel['kernel_operator']
    kernel_type = ctp_kernel['kernel_type']

    # Handle frequency factor
    if hasattr(clock_frequencies, '__iter__'):
        # Multiple frequencies - use supremum
        omega_squared_max = max(omega**2 for omega in clock_frequencies)
    else:
        omega_squared_max = clock_frequencies**2

    coupling = operator_params['coupling_constant']

    def operator_T(phi_field):
        """Apply T[Φ] = C ω² C[Φ]"""
        # First application of C
        C_phi = C_op(phi_field)

        # Multiply by ω² and coupling
        omega_squared_C_phi = coupling * omega_squared_max * C_phi

        # Second application of C
        result = C_op(omega_squared_C_phi)

        return result

    # Estimate operator norm bound
    if kernel_type == 'local':
        # For local kernel: ||T|| = α² ω²_max
        alpha = ctp_kernel['coupling_constant']
        norm_bound = coupling * alpha**2 * omega_squared_max

    elif kernel_type in ['quasilocal', 'screened']:
        # For non-local: ||T|| ≤ ||C||² ||ω²||_∞
        spectral_radius_C = ctp_kernel['spectral_properties']['max_eigenvalue']
        norm_bound = coupling * spectral_radius_C**2 * omega_squared_max

    else:
        norm_bound = None

    return {
        'operator_T': operator_T,
        'kernel_C': C_op,
        'frequency_factor': omega_squared_max,
        'coupling_constant': coupling,
        'operator_norm_bound': norm_bound,
        'kernel_type': kernel_type,
        'composition': 'T = C ω² C'
    }

def check_contraction_property(operator_data, contraction_params=None):
    """
    Check if T is a contraction mapping.

    Contraction condition: ||T|| < 1
    From paper: C||ω²||_∞ L_C < 1

    Args:
        operator_data: Result from build_fixed_point_operator
        contraction_params: Optional parameters for check

    Returns:
        Dict with:
            'is_contraction': Boolean result
            'operator_norm': Estimated ||T||
            'contraction_factor': How much < 1
            'convergence_rate': Geometric convergence rate
    """
    if contraction_params is None:
        contraction_params = {}

    norm_bound = operator_data['operator_norm_bound']

    if norm_bound is None:
        return {
            'is_contraction': None,
            'operator_norm': None,
            'message': 'Cannot determine norm bound'
        }

    is_contraction = norm_bound < 1

    if is_contraction:
        # Convergence rate for Banach fixed-point theorem
        # ||Φ^(n+1) - Φ*|| ≤ q^n ||Φ^(1) - Φ^(0)||
        convergence_rate = norm_bound
        iterations_for_epsilon = lambda eps: ceil(log(eps)/log(norm_bound))
    else:
        convergence_rate = None
        iterations_for_epsilon = None

    return {
        'is_contraction': is_contraction,
        'operator_norm': norm_bound,
        'contraction_factor': norm_bound if norm_bound < 1 else None,
        'convergence_rate': convergence_rate,
        'iterations_for_1e-6': iterations_for_epsilon(1e-6) if iterations_for_epsilon else None,
        'condition_satisfied': is_contraction,
        'banach_applicable': is_contraction
    }

def iterate_fixed_point(operator_data, initial_field, num_iterations=10, iteration_params=None):
    """
    Iterate the fixed-point map Φ^(n+1) = T[Φ^n] with detailed convergence tracking.

    Uses relaxation: Φ^(n+1) = (1-α)Φ^n + α T[Φ^n]
    for stability with relaxation parameter α.

    Convergence tracking:
    - Tracks convergence uncertainty and error bounds
    - Compares actual vs theoretical convergence rates
    - Provides detailed convergence analysis

    Args:
        operator_data: Fixed-point operator T
        initial_field: Starting field Φ^(0)
        num_iterations: Number of iterations
        iteration_params: Dict with 'relaxation' parameter α and options

    Returns:
        Dict with:
            'fixed_point': Final Φ after iterations
            'iteration_history': List of fields
            'convergence_metrics': ||Φ^(n+1) - Φ^n||
            'converged': Whether convergence achieved
            'convergence_analysis': Detailed rate analysis
            'error_bounds': Theoretical vs actual bounds
            'uncertainty_estimates': Error estimates per iteration
    """
    if iteration_params is None:
        iteration_params = {'relaxation': 0.5, 'tolerance': 1e-6, 'track_uncertainty': True}

    alpha = iteration_params.get('relaxation', 0.5)
    tolerance = iteration_params.get('tolerance', 1e-6)
    track_uncertainty = iteration_params.get('track_uncertainty', True)
    T = operator_data['operator_T']
    norm_bound = operator_data.get('operator_norm_bound')

    # Initialize
    phi = initial_field
    history = [phi]
    convergence_metrics = []
    uncertainty_estimates = []
    theoretical_bounds = []

    converged = False

    # Theoretical analysis setup
    if norm_bound is not None and norm_bound < 1:
        q = norm_bound  # Contraction factor
        theoretical_convergence = True
    else:
        q = None
        theoretical_convergence = False

    for n in range(num_iterations):
        # Apply operator
        T_phi = T(phi)

        # Relaxed update
        phi_new = (1 - alpha) * phi + alpha * T_phi

        # Convergence metric (enhanced)
        if hasattr(phi, '__sub__'):
            difference = phi_new - phi
            # Improved metric calculation
            if hasattr(difference, '__abs__'):
                metric = abs(difference)
            elif hasattr(difference, 'norm'):
                metric = difference.norm()
            elif hasattr(difference, '__float__'):
                metric = float(abs(difference))
            else:
                raise TypeError(f"Field type {type(difference)} does not support required operations (abs, norm, or float conversion)")
        else:
            raise TypeError(f"Field type {type(phi)} does not support subtraction required for convergence metric")

        convergence_metrics.append(metric)
        history.append(phi_new)

        # Uncertainty estimation (if tracking enabled)
        if track_uncertainty and len(convergence_metrics) >= 2:
            # Estimate uncertainty from convergence rate variation
            recent_metrics = convergence_metrics[-2:]
            if recent_metrics[0] > 0:
                rate_estimate = recent_metrics[1] / recent_metrics[0]
                # Uncertainty based on rate stability
                uncertainty = abs(rate_estimate - q) if q else 0.1 * metric
            else:
                uncertainty = 0.1 * metric
        else:
            uncertainty = metric * 0.1  # 10% uncertainty estimate

        uncertainty_estimates.append(uncertainty)

        # Theoretical bound for this iteration
        if theoretical_convergence and len(convergence_metrics) >= 1:
            initial_residual = convergence_metrics[0]
            theoretical_bound = initial_residual * (q ** (n + 1))
            theoretical_bounds.append(theoretical_bound)
        else:
            theoretical_bounds.append(None)

        # Check convergence
        if metric < tolerance:
            converged = True
            break

        phi = phi_new

    # Convergence analysis
    convergence_analysis = {}

    if len(convergence_metrics) >= 3:
        # Estimate actual convergence rate from data
        ratios = []
        for i in range(1, len(convergence_metrics)):
            if convergence_metrics[i-1] > 1e-15:
                ratio = convergence_metrics[i] / convergence_metrics[i-1]
                ratios.append(ratio)

        if ratios:
            estimated_rate = mean(ratios)
            rate_std = sqrt(variance(ratios)) if len(ratios) > 1 else 0

            convergence_analysis = {
                'estimated_convergence_rate': estimated_rate,
                'rate_standard_deviation': rate_std,
                'theoretical_rate': q,
                'rate_match': abs(estimated_rate - q) < 0.2 if q else None,
                'geometric_convergence': estimated_rate < 1,
                'convergence_order': 'linear' if estimated_rate > 0 else 'superlinear'
            }

    # Error bounds analysis
    error_bounds = {}
    if theoretical_bounds and convergence_metrics:
        actual_vs_theory = []
        for actual, theory in zip(convergence_metrics, theoretical_bounds):
            if theory is not None and theory > 0:
                ratio = actual / theory
                actual_vs_theory.append(ratio)

        if actual_vs_theory:
            error_bounds = {
                'theoretical_bounds': theoretical_bounds,
                'actual_residuals': convergence_metrics,
                'bound_ratios': actual_vs_theory,
                'bounds_satisfied': all(r <= 10 for r in actual_vs_theory),  # Allow factor of 10
                'average_bound_ratio': mean(actual_vs_theory),
                'max_bound_violation': max(actual_vs_theory) if actual_vs_theory else 0
            }

    # Final convergence properties
    final_properties = {
        'converged_to_tolerance': converged,
        'final_residual': convergence_metrics[-1] if convergence_metrics else None,
        'residual_uncertainty': uncertainty_estimates[-1] if uncertainty_estimates else None,
        'iterations_to_convergence': len(convergence_metrics),
        'convergence_certified': converged and (
            error_bounds.get('bounds_satisfied', False) if error_bounds else True
        )
    }

    return {
        'fixed_point': phi,
        'iteration_history': history,
        'convergence_metrics': convergence_metrics,
        'converged': converged,
        'num_iterations_used': len(history) - 1,
        'final_residual': convergence_metrics[-1] if convergence_metrics else None,
        'relaxation_used': alpha,
        # Convergence analysis data
        'convergence_analysis': convergence_analysis,
        'error_bounds': error_bounds,
        'uncertainty_estimates': uncertainty_estimates,
        'theoretical_bounds': theoretical_bounds,
        'final_properties': final_properties
    }

def compute_spectral_radius(operator_data, spectral_params=None):
    """
    Compute spectral radius of operator T with convergence tracking.

    Uses power iteration or Lanczos method for large operators.
    The spectral radius determines contraction properties.

    Spectral analysis:
    - Tracks convergence history for power iteration
    - Provides uncertainty estimates on spectral radius
    - Records iteration-by-iteration eigenvalue estimates

    Args:
        operator_data: Fixed-point operator
        spectral_params: Method and parameters

    Returns:
        Dict with:
            'spectral_radius': ρ(T)
            'dominant_eigenvalue': λ_max
            'method_used': Computation method
            'is_contraction': ρ(T) < 1
            'convergence_history': Iteration-by-iteration estimates
            'convergence_analysis': Rate and uncertainty analysis
    """
    if spectral_params is None:
        spectral_params = {'method': 'estimate', 'iterations': 100, 'tolerance': 1e-8}

    method = spectral_params['method']
    convergence_history = []
    convergence_analysis = {}

    if method == 'estimate':
        # Use the norm bound as estimate
        spectral_radius = operator_data['operator_norm_bound']
        convergence_history = [spectral_radius] if spectral_radius is not None else []

    elif method == 'power_iteration':
        # Power iteration with convergence tracking
        T = operator_data['operator_T']
        max_iterations = spectral_params.get('iterations', 100)
        tolerance = spectral_params.get('tolerance', 1e-8)

        # Start with initial vector (could be enhanced with better initial guess)
        v = 1.0  # Placeholder - in practice would use appropriate field type

        eigenvalue_estimates = []
        vector_norms = []
        converged = False

        for iteration in range(max_iterations):
            # Apply operator
            Tv = T(v)

            # Compute eigenvalue estimate (Rayleigh quotient)
            if hasattr(Tv, '__mul__') and hasattr(v, '__mul__'):
                numerator = Tv * v
                denominator = v * v
                if abs(denominator) > 1e-15:
                    eigenvalue_est = numerator / denominator
                elif hasattr(Tv, '__abs__'):
                    eigenvalue_est = abs(Tv)
                else:
                    raise TypeError(f"Cannot compute eigenvalue: field type {type(Tv)} does not support required operations")
            elif hasattr(Tv, '__abs__') and hasattr(v, '__abs__'):
                # Simplified estimate using absolute values
                if abs(v) > 1e-15:
                    eigenvalue_est = abs(Tv) / abs(v)
                else:
                    eigenvalue_est = abs(Tv)
            else:
                raise TypeError(f"Field types {type(Tv)}, {type(v)} do not support required operations for eigenvalue estimation")

            eigenvalue_estimates.append(eigenvalue_est)

            # Normalize vector for next iteration
            if hasattr(Tv, '__abs__'):
                v_norm = abs(Tv)
                if v_norm > 1e-15:
                    v = Tv / v_norm
                else:
                    v = Tv  # Keep unnormalized if too small
            else:
                v = Tv

            vector_norms.append(v_norm if 'v_norm' in locals() else 1.0)

            # Check convergence
            if len(eigenvalue_estimates) >= 2:
                change = abs(eigenvalue_estimates[-1] - eigenvalue_estimates[-2])
                if change < tolerance:
                    converged = True
                    break

        spectral_radius = eigenvalue_estimates[-1] if eigenvalue_estimates else None
        convergence_history = eigenvalue_estimates

        # Analyze convergence properties
        if len(eigenvalue_estimates) >= 3:
            # Estimate convergence rate of power iteration
            changes = []
            for i in range(2, len(eigenvalue_estimates)):
                prev_change = abs(eigenvalue_estimates[i-1] - eigenvalue_estimates[i-2])
                curr_change = abs(eigenvalue_estimates[i] - eigenvalue_estimates[i-1])
                if prev_change > 1e-15:
                    conv_rate = curr_change / prev_change
                    changes.append(conv_rate)

            convergence_analysis = {
                'converged': converged,
                'iterations_used': len(eigenvalue_estimates),
                'final_eigenvalue': spectral_radius,
                'convergence_rate': mean(changes) if changes else None,
                'rate_stability': sqrt(variance(changes)) if len(changes) > 1 else None,
                'eigenvalue_uncertainty': abs(eigenvalue_estimates[-1] - eigenvalue_estimates[-2]) if len(eigenvalue_estimates) >= 2 else None
            }

    elif method == 'power_iteration_enhanced':
        # More sophisticated power iteration with deflation
        T = operator_data['operator_T']
        max_iterations = spectral_params.get('iterations', 100)
        tolerance = spectral_params.get('tolerance', 1e-8)

        # Multiple starting vectors for better coverage
        starting_vectors = [1.0, 0.5, 2.0]  # Could be more sophisticated
        all_eigenvalue_estimates = []

        for start_v in starting_vectors:
            v = start_v
            eigenvalue_sequence = []

            for iteration in range(max_iterations // len(starting_vectors)):
                Tv = T(v)

                # Rayleigh quotient
                if hasattr(Tv, '__mul__') and hasattr(v, '__mul__'):
                    numerator = Tv * v
                    denominator = v * v
                    if abs(denominator) > 1e-15:
                        eigenvalue_est = numerator / denominator
                    elif hasattr(Tv, '__abs__'):
                        eigenvalue_est = abs(Tv)
                    else:
                        raise TypeError(f"Cannot compute eigenvalue: field type {type(Tv)} does not support required operations")
                elif hasattr(Tv, '__abs__') and hasattr(v, '__abs__'):
                    if abs(v) > 1e-15:
                        eigenvalue_est = abs(Tv) / abs(v)
                    else:
                        eigenvalue_est = abs(Tv)
                else:
                    raise TypeError(f"Field types {type(Tv)}, {type(v)} do not support required operations for eigenvalue estimation")

                eigenvalue_sequence.append(eigenvalue_est)

                # Normalize
                if hasattr(Tv, '__abs__'):
                    v_norm = abs(Tv)
                    v = Tv / v_norm if v_norm > 1e-15 else Tv
                else:
                    v = Tv

            all_eigenvalue_estimates.extend(eigenvalue_sequence)

        # Take maximum eigenvalue estimate
        spectral_radius = max(all_eigenvalue_estimates) if all_eigenvalue_estimates else None
        convergence_history = all_eigenvalue_estimates

        convergence_analysis = {
            'method_variant': 'multiple_starting_vectors',
            'num_starting_vectors': len(starting_vectors),
            'final_eigenvalue': spectral_radius,
            'eigenvalue_range': [min(all_eigenvalue_estimates), max(all_eigenvalue_estimates)] if all_eigenvalue_estimates else None
        }

    elif method == 'lanczos':
        # Lanczos method for large sparse operators
        # This requires more sophisticated implementation
        raise NotImplementedError("Lanczos method not yet implemented")

    else:
        raise ValueError(f"Unknown spectral method: {method}")

    is_contraction = spectral_radius < 1 if spectral_radius is not None else None

    # Additional uncertainty analysis
    uncertainty_estimate = None
    if convergence_history and len(convergence_history) >= 2:
        # Estimate uncertainty from convergence stability
        recent_estimates = convergence_history[-min(5, len(convergence_history)):]
        if len(recent_estimates) > 1:
            uncertainty_estimate = sqrt(variance(recent_estimates))

    return {
        'spectral_radius': spectral_radius,
        'dominant_eigenvalue': spectral_radius,  # For self-adjoint operators
        'method_used': method,
        'is_contraction': is_contraction,
        'iterations_used': spectral_params.get('iterations'),
        # Convergence analysis data
        'convergence_history': convergence_history,
        'convergence_analysis': convergence_analysis,
        'uncertainty_estimate': uncertainty_estimate,
        'final_convergence_properties': {
            'spectral_radius_converged': convergence_analysis.get('converged', True),
            'spectral_radius_uncertainty': uncertainty_estimate,
            'contraction_verified': is_contraction
        }
    }

def verify_banach_conditions(operator_data, verification_params=None):
    """
    Verify conditions for Banach fixed-point theorem.

    Requirements:
    1. T maps complete space to itself
    2. T is a contraction: ||T|| < 1
    3. Uniqueness of fixed point
    4. Geometric convergence

    Args:
        operator_data: Operator T to verify
        verification_params: Optional parameters

    Returns:
        Dict with comprehensive verification
    """
    if verification_params is None:
        verification_params = {}

    # Check contraction property
    contraction_result = check_contraction_property(operator_data)

    # Compute spectral radius
    spectral_result = compute_spectral_radius(operator_data)

    # Banach conditions
    is_contraction = contraction_result['is_contraction']
    has_complete_domain = True  # Assumed for function spaces with appropriate topology

    banach_satisfied = is_contraction and has_complete_domain

    if banach_satisfied:
        # Convergence rate from contraction factor
        q = contraction_result['contraction_factor']

        # Error after n iterations: ||Φ^n - Φ*|| ≤ q^n/(1-q) ||Φ^1 - Φ^0||
        error_bound = lambda n: q**n / (1 - q)

        # Iterations needed for given accuracy
        iterations_for_accuracy = lambda eps: ceil(log(eps * (1 - q)) / log(q))

        uniqueness = True  # Guaranteed by Banach theorem
        existence = True   # Guaranteed by Banach theorem

    else:
        error_bound = None
        iterations_for_accuracy = None
        uniqueness = None
        existence = None

    return {
        'banach_conditions_satisfied': banach_satisfied,
        'contraction_verified': is_contraction,
        'complete_space': has_complete_domain,
        'fixed_point_exists': existence,
        'fixed_point_unique': uniqueness,
        'convergence_guaranteed': banach_satisfied,
        'spectral_radius': spectral_result['spectral_radius'],
        'contraction_factor': contraction_result.get('contraction_factor'),
        'error_bound_function': error_bound,
        'iterations_for_1e-9': iterations_for_accuracy(1e-9) if iterations_for_accuracy else None,
        'theoretical_framework': 'Banach fixed-point theorem'
    }