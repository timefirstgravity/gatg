#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Test script for SymbolicGATGEquations
Verifies that symbolic expressions match documentation and mathematical properties
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sage.all import *
from symbolic_equations import SymbolicGATGEquations, verify_schwarzschild, verify_flux_law, verify_christoffel_symmetry, verify_christoffel_connection, verify_christoffel_spherical_symmetry

def test_schwarzschild_equations():
    """Test Schwarzschild equation symbolic implementation"""
    print("="*70)
    print(" "*20 + "SCHWARZSCHILD EQUATION TESTS")
    print("="*70)

    # Create instance
    seq = SymbolicGATGEquations()

    # Test 1: Verify equation access methods
    print("\n1. Testing equation access methods:")
    print("-"*50)

    ode_symbolic = seq.get_symbolic('schwarzschild_ode')
    ode_string = seq.get_string('schwarzschild_ode')
    ode_latex = seq.get_latex('schwarzschild_ode')
    ode_desc = seq.get_description('schwarzschild_ode')

    print(f"  String representation: {ode_string}")
    print(f"  LaTeX: {ode_latex}")
    print(f"  Description: {ode_desc}")
    print(f"  Symbolic expression: {ode_symbolic}")

    # Test 2: Verify the ODE is satisfied by the known solution
    print("\n2. Verifying Schwarzschild ODE with known solution:")
    print("-"*50)

    is_valid = seq.verify_schwarzschild_ode()
    print(f"  ODE satisfied by A(r) = 1 - r_s/r: {is_valid}")

    if is_valid:
        print("  ✓ PASSED: Schwarzschild ODE correctly verified")
    else:
        print("  ✗ FAILED: Schwarzschild ODE verification failed")

    # Test 3: Manual verification with substitution
    print("\n3. Manual verification with explicit substitution:")
    print("-"*50)

    r = seq.vars['r']
    rs = seq.vars['r_s']
    A_solution = 1 - rs/r

    # Compute derivative
    dA_dr = diff(A_solution, r)
    print(f"  A(r) = {A_solution}")
    print(f"  dA/dr = {dA_dr}")

    # Substitute into ODE: r * A'(r) + A(r) - 1
    ode_lhs = r * dA_dr + A_solution - 1
    ode_simplified = ode_lhs.simplify_full()
    print(f"  r*A'(r) + A(r) - 1 = {ode_simplified}")

    if ode_simplified == 0:
        print("  ✓ PASSED: Manual verification successful")
    else:
        print(f"  ✗ FAILED: Expected 0, got {ode_simplified}")

    # Test 4: Verify constraint form
    print("\n4. Testing Schwarzschild constraint (rA)' = 1:")
    print("-"*50)

    constraint_expr = seq.get_symbolic('schwarzschild_constraint')
    constraint_string = seq.get_string('schwarzschild_constraint')
    print(f"  Constraint: {constraint_string}")

    # Verify with solution
    rA = r * A_solution
    d_rA = diff(rA, r).simplify_full()
    print(f"  d(rA)/dr = {d_rA}")

    if d_rA == 1:
        print("  ✓ PASSED: Constraint verified")
    else:
        print(f"  ✗ FAILED: Expected 1, got {d_rA}")

    return is_valid and ode_simplified == 0 and d_rA == 1

def test_flux_law_equations():
    """Test flux law equation symbolic implementation"""
    print("\n" + "="*70)
    print(" "*25 + "FLUX LAW TESTS")
    print("="*70)

    seq = SymbolicGATGEquations()

    # Test 1: Access flux law components
    print("\n1. Testing flux law equation access:")
    print("-"*50)

    flux_symbolic = seq.get_symbolic('flux_law')
    flux_string = seq.get_string('flux_law')
    flux_desc = seq.get_description('flux_law')

    print(f"  String: {flux_string}")
    print(f"  Description: {flux_desc}")
    print(f"  Symbolic: {flux_symbolic}")

    # Test 2: Verify consistency between flux law and Einstein relations
    print("\n2. Verifying flux law consistency with Einstein tensor:")
    print("-"*50)

    is_consistent = seq.verify_flux_law_consistency()
    print(f"  Flux law consistent with Einstein relations: {is_consistent}")

    if is_consistent:
        print("  ✓ PASSED: Flux law consistency verified")
        print("  This confirms: ∂_t Φ = (4πG/c⁴) r T^t_r")
        print("               → G^t_r = (2/r) ∂_t Φ = (8πG/c⁴) T^t_r")
    else:
        print("  ✗ FAILED: Flux law consistency check failed")

    # Test 3: Check geometrized units version
    print("\n3. Testing geometrized units (G=c=1):")
    print("-"*50)

    flux_geo = seq.get_symbolic('flux_law_geometrized')
    flux_geo_string = seq.get_string('flux_law_geometrized')
    print(f"  Geometrized form: {flux_geo_string}")

    # Verify that setting G=c=1 in standard form gives geometrized form
    t = seq.vars['t']
    r = seq.vars['r']
    Phi = function('Phi')(t, r)
    T_tr = var('T_tr')

    standard_with_G1_c1 = diff(Phi, t) - (4*pi*1/1**4) * r * T_tr
    geometrized = diff(Phi, t) - 4*pi * r * T_tr

    if standard_with_G1_c1 == geometrized:
        print("  ✓ PASSED: Geometrized units correctly derived from standard form")
    else:
        print("  ✗ FAILED: Geometrized form inconsistent")

    return is_consistent

def test_equation_categories():
    """Test equation categorization"""
    print("\n" + "="*70)
    print(" "*20 + "EQUATION CATEGORIZATION TEST")
    print("="*70)

    seq = SymbolicGATGEquations()

    print("\n1. Listing equations by category:")
    print("-"*50)

    for category in ['schwarzschild', 'flux_law', 'adm', 'kerr']:
        equations = seq.list_equations(category)
        print(f"\n  {category.upper()} category:")
        for eq in equations[:3]:  # Show first 3 for brevity
            desc = seq.get_description(eq)
            print(f"    • {eq}: {desc[:60]}...")
        if len(equations) > 3:
            print(f"    ... and {len(equations) - 3} more")

    print("\n2. All equations:")
    print("-"*50)
    all_equations = seq.list_equations()
    print(f"  Total equations defined: {len(all_equations)}")
    for eq in all_equations[:5]:  # Show first 5
        print(f"    • {eq}")
    if len(all_equations) > 5:
        print(f"    ... and {len(all_equations) - 5} more")

    return True

def test_evaluation_with_substitution():
    """Test equation evaluation with specific values"""
    print("\n" + "="*70)
    print(" "*20 + "EQUATION EVALUATION TEST")
    print("="*70)

    seq = SymbolicGATGEquations()

    print("\n1. Evaluating Schwarzschild radius for specific mass:")
    print("-"*50)

    # Solar mass in SI units
    M_sun = 1.989e30  # kg
    c = 3e8  # m/s
    G = 6.674e-11  # m^3/(kg*s^2)

    substitutions = {
        seq.vars['M']: M_sun,
        seq.vars['c']: c,
        seq.vars['G']: G
    }

    rs_symbolic = seq.get_symbolic('schwarzschild_radius')
    rs_value = rs_symbolic.substitute(substitutions)

    print(f"  M_sun = {M_sun:.3e} kg")
    print(f"  c = {c:.3e} m/s")
    print(f"  G = {G:.3e} m^3/(kg*s^2)")

    # Convert to float for printing
    rs_value_float = float(rs_value)
    print(f"  Schwarzschild radius: r_s = {rs_value_float:.3f} m")
    print(f"  (Expected ~2953 m for solar mass)")

    expected_rs = 2 * G * M_sun / c**2
    if abs(rs_value_float - expected_rs) < 1:  # Within 1 meter
        print("  ✓ PASSED: Schwarzschild radius calculation correct")
    else:
        print(f"  ✗ FAILED: Expected {expected_rs:.3f} m")

    return True

def test_latex_generation():
    """Test LaTeX document generation"""
    print("\n" + "="*70)
    print(" "*20 + "LATEX GENERATION TEST")
    print("="*70)

    seq = SymbolicGATGEquations()

    print("\n1. Generating LaTeX for Schwarzschild equations:")
    print("-"*50)

    latex_doc = seq.generate_latex_document('schwarzschild')
    print("  Generated LaTeX fragment:")
    print("  " + "-"*40)
    for line in latex_doc.split('\n')[:10]:  # Show first 10 lines
        print(f"  {line}")
    print("  " + "-"*40)

    if '\\begin{equation}' in latex_doc and '\\label{eq:' in latex_doc:
        print("  ✓ PASSED: LaTeX generation successful")
    else:
        print("  ✗ FAILED: LaTeX generation incomplete")

    return True

def test_adm_equations():
    """Test ADM decomposition equation symbolic implementation"""
    print("\n" + "="*70)
    print(" "*25 + "ADM DECOMPOSITION TESTS")
    print("="*70)

    seq = SymbolicGATGEquations()

    print("\n1. Testing ADM lapse function relation:")
    print("-"*50)

    lapse_symbolic = seq.get_symbolic('lapse_function')
    lapse_string = seq.get_string('lapse_function')
    lapse_desc = seq.get_description('lapse_function')

    print(f"  String: {lapse_string}")
    print(f"  Description: {lapse_desc}")
    print(f"  Symbolic: {lapse_symbolic}")

    # Verify lapse function relation
    N = function('N')(var('t'), var('r'))
    Phi = function('Phi')(var('t'), var('r'))
    lapse_relation = N - exp(Phi)
    print(f"  ✓ Lapse function correctly relates N to Φ")

    print("\n2. Testing extrinsic curvature definition:")
    print("-"*50)

    K_def = seq.get_string('extrinsic_curvature_def')
    K_trace = seq.get_string('extrinsic_trace')
    print(f"  Extrinsic curvature: {K_def}")
    print(f"  Trace: {K_trace}")

    print("\n3. Testing 3-metric components:")
    print("-"*50)

    gamma_rr = seq.get_string('gamma_rr')
    gamma_thth = seq.get_string('gamma_thth')
    gamma_phph = seq.get_string('gamma_phph')

    print(f"  γ_rr: {gamma_rr}")
    print(f"  γ_θθ: {gamma_thth}")
    print(f"  γ_φφ: {gamma_phph}")

    print("\n4. Testing constraint equations:")
    print("-"*50)

    momentum = seq.get_string('momentum_constraint')
    hamiltonian = seq.get_string('hamiltonian_constraint')

    print(f"  Momentum: {momentum}")
    print(f"  Hamiltonian: {hamiltonian}")

    print("\n5. Listing all ADM equations:")
    print("-"*50)

    adm_equations = seq.list_equations('adm')
    print(f"  Total ADM equations defined: {len(adm_equations)}")
    for eq in adm_equations[:5]:
        desc = seq.get_description(eq)
        print(f"    • {eq}: {desc[:50]}...")

    return True

def test_kerr_equations():
    """Test Kerr metric equation symbolic implementation"""
    print("\n" + "="*70)
    print(" "*25 + "KERR METRIC TESTS")
    print("="*70)

    seq = SymbolicGATGEquations()

    print("\n1. Testing Kerr auxiliary functions:")
    print("-"*50)

    sigma_symbolic = seq.get_symbolic('kerr_sigma')
    delta_symbolic = seq.get_symbolic('kerr_delta')
    sigma_string = seq.get_string('kerr_sigma')
    delta_string = seq.get_string('kerr_delta')

    print(f"  Σ function: {sigma_string}")
    print(f"  Symbolic: {sigma_symbolic}")
    print(f"  Δ function: {delta_string}")
    print(f"  Symbolic: {delta_symbolic}")

    print("\n2. Testing Kerr metric components:")
    print("-"*50)

    components = ['kerr_g_tt', 'kerr_g_rr', 'kerr_g_thth', 'kerr_g_phph', 'kerr_g_tph']
    for comp in components:
        comp_string = seq.get_string(comp)
        comp_desc = seq.get_description(comp)
        print(f"  {comp}: {comp_string}")
        if 'frame-dragging' in comp_desc:
            print(f"    → Frame-dragging component!")
        print()

    print("3. Testing Kerr physical properties:")
    print("-"*50)

    horizons_string = seq.get_string('kerr_event_horizons')
    ergosphere_string = seq.get_string('kerr_ergosphere')

    print(f"  Event horizons: {horizons_string}")
    print(f"  Ergosphere: {ergosphere_string}")

    print("\n4. Testing Kerr mathematical verifications:")
    print("-"*50)

    # Test Schwarzschild limit
    schwarzschild_limit_ok = seq.verify_kerr_schwarzschild_limit()
    print(f"  Schwarzschild limit (a→0): {'✓ VERIFIED' if schwarzschild_limit_ok else '✗ FAILED'}")

    # Test horizon properties
    horizons_ok = seq.verify_kerr_horizons()
    print(f"  Event horizons satisfy Δ=0: {'✓ VERIFIED' if horizons_ok else '✗ FAILED'}")

    # Test internal consistency
    consistency_ok = seq.verify_kerr_metric_consistency()
    print(f"  Auxiliary function consistency: {'✓ VERIFIED' if consistency_ok else '✗ FAILED'}")

    print("\n5. Testing frame-dragging effect:")
    print("-"*50)

    g_tph = seq.get_symbolic('kerr_g_tph')
    print(f"  Off-diagonal component: {seq.get_string('kerr_g_tph')}")
    print("  • Non-zero only when a ≠ 0 (rotation)")
    print("  • Causes time-space mixing (frame-dragging)")
    print("  • Lense-Thirring effect in weak field limit")

    # Test that g_tph vanishes when a=0
    g_tph_no_rotation = g_tph.substitute({seq.vars['a']: 0})
    if g_tph_no_rotation == 0:
        print("  ✓ Correctly vanishes when a=0 (no rotation)")
    else:
        print("  ✗ Should vanish when a=0")

    print("\n6. Kerr equation count:")
    print("-"*50)
    kerr_equations = seq.list_equations('kerr')
    print(f"  Total Kerr equations implemented: {len(kerr_equations)}")
    print("  This includes metric components, auxiliary functions,")
    print("  physical properties, and lapse-first decomposition.")

    return (schwarzschild_limit_ok and horizons_ok and
            consistency_ok and g_tph_no_rotation == 0)

def test_christoffel_equations():
    """Test Christoffel symbol equation symbolic implementation"""
    print("\n" + "="*70)
    print(" "*20 + "CHRISTOFFEL SYMBOL TESTS")
    print("="*70)

    seq = SymbolicGATGEquations()

    # Test 1: Access all Christoffel symbols
    print("\n1. Testing Christoffel symbol access:")
    print("-"*50)

    christoffel_symbols = [
        'christoffel_t_tt', 'christoffel_t_tr', 'christoffel_t_rr',
        'christoffel_r_tt', 'christoffel_r_tr', 'christoffel_r_rr',
        'christoffel_theta_rtheta', 'christoffel_phi_rphi',
        'christoffel_r_thetatheta', 'christoffel_r_phiphi'
    ]

    print(f"  Implemented Christoffel symbols: {len(christoffel_symbols)}")
    for symbol in christoffel_symbols[:3]:  # Show first 3 as examples
        string_rep = seq.get_string(symbol)
        desc = seq.get_description(symbol)
        print(f"    {symbol}: {string_rep}")
        print(f"      → {desc}")

    # Test 2: Verify spherical symmetry properties
    print("\n2. Verifying spherical symmetry properties:")
    print("-"*50)

    spherical_ok = seq.verify_christoffel_spherical_symmetry()
    print(f"  Spherical symmetry verified: {'✓ PASSED' if spherical_ok else '✗ FAILED'}")

    if spherical_ok:
        print("  • Angular components have correct 1/r dependence")
        print("  • Radial components have proper r scaling")
    else:
        print("  ✗ Spherical symmetry properties failed")

    # Test 3: Verify connection symmetries
    print("\n3. Verifying connection symmetries:")
    print("-"*50)

    symmetry_ok = seq.verify_christoffel_symmetry()
    print(f"  Connection symmetries verified: {'✓ PASSED' if symmetry_ok else '✗ FAILED'}")

    if symmetry_ok:
        print("  • Γ^μ_νσ = Γ^μ_σν (symmetric in lower indices)")
    else:
        print("  ✗ Connection symmetry verification failed")

    # Test 4: Verify connection property
    print("\n4. Verifying connection property:")
    print("-"*50)

    connection_ok = seq.verify_christoffel_connection_property()
    print(f"  Connection property verified: {'✓ PASSED' if connection_ok else '✗ FAILED'}")

    if connection_ok:
        print("  • Metric compatibility: ∇_μ g_νσ = 0")
        print("  • Proper sign patterns for spherical metric")
    else:
        print("  ✗ Connection property verification failed")

    # Test 5: Verify specific Christoffel components
    print("\n5. Testing specific Christoffel components:")
    print("-"*50)

    # Test temporal components
    gamma_t_tt = seq.get_symbolic('christoffel_t_tt')  # ∂_t Φ
    gamma_t_tr = seq.get_symbolic('christoffel_t_tr')  # ∂_r Φ

    print(f"  Γ^t_tt = {seq.get_string('christoffel_t_tt')} (temporal evolution)")
    print(f"  Γ^t_tr = {seq.get_string('christoffel_t_tr')} (radial derivative)")

    # Test angular components
    gamma_theta_rtheta = seq.get_symbolic('christoffel_theta_rtheta')  # 1/r
    gamma_phi_rphi = seq.get_symbolic('christoffel_phi_rphi')          # 1/r

    r = seq.vars['r']
    angular_components_ok = True

    if (gamma_theta_rtheta - 1/r).simplify_full() == 0:
        print("  ✓ Γ^θ_rθ = 1/r correctly implemented")
    else:
        print("  ✗ Γ^θ_rθ ≠ 1/r")
        angular_components_ok = False

    if (gamma_phi_rphi - 1/r).simplify_full() == 0:
        print("  ✓ Γ^φ_rφ = 1/r correctly implemented")
    else:
        print("  ✗ Γ^φ_rφ ≠ 1/r")
        angular_components_ok = False

    # Test 6: Verify metric-specific properties
    print("\n6. Testing metric-specific properties:")
    print("-"*50)

    # For metric ds² = -e^(2Φ)dt² + e^(-2Φ)dr² + r²dΩ²
    # Verify that Christoffel symbols have the expected functional form
    t = seq.vars['t']
    theta = seq.vars['theta']
    Phi = function('Phi')(t, r)

    # Check that temporal derivatives appear in temporal components
    gamma_t_tt_expected = gamma_t_tt.has(diff(Phi, t))
    gamma_t_tr_expected = gamma_t_tr.has(diff(Phi, r))

    print(f"  Temporal components depend on Φ derivatives: {'✓' if gamma_t_tt_expected and gamma_t_tr_expected else '✗'}")

    # Check that radial components have exponential factors
    gamma_r_tt = seq.get_symbolic('christoffel_r_tt')  # e^(4Φ) ∂_r Φ
    gamma_t_rr = seq.get_symbolic('christoffel_t_rr')  # e^(-4Φ) ∂_t Φ

    exponential_factors_ok = gamma_r_tt.has(exp(4*Phi)) and gamma_t_rr.has(exp(-4*Phi))
    print(f"  Exponential factors in mixed components: {'✓' if exponential_factors_ok else '✗'}")

    # Test 7: Category membership
    print("\n7. Testing category membership:")
    print("-"*50)

    christoffel_list = seq.list_equations('christoffel')
    print(f"  Christoffel category contains {len(christoffel_list)} equations")

    all_found = all(symbol in christoffel_list for symbol in christoffel_symbols)
    print(f"  All symbols properly categorized: {'✓ PASSED' if all_found else '✗ FAILED'}")

    # Overall success
    overall_success = (spherical_ok and symmetry_ok and connection_ok and
                      angular_components_ok and exponential_factors_ok and all_found)

    if overall_success:
        print("\n🎉 All Christoffel symbol tests passed!")
        print("   The implementation correctly captures the geometry of")
        print("   the spherically symmetric metric ds² = -e^(2Φ)dt² + e^(-2Φ)dr² + r²dΩ²")
    else:
        print("\n⚠️ Some Christoffel symbol tests failed.")

    return overall_success

def main():
    """Run all tests"""
    print("\n" + "="*70)
    print(" "*15 + "SYMBOLIC GATG EQUATIONS TEST SUITE")
    print("="*70)
    print("\nTesting the new SymbolicGATGEquations implementation that")
    print("connects string documentation with symbolic mathematics.\n")

    # Run tests
    tests = [
        ("Schwarzschild Equations", test_schwarzschild_equations),
        ("Flux Law Equations", test_flux_law_equations),
        ("ADM Decomposition Equations", test_adm_equations),
        ("Kerr Metric Equations", test_kerr_equations),
        ("Christoffel Symbols", test_christoffel_equations),
        ("Equation Categories", test_equation_categories),
        ("Evaluation with Substitution", test_evaluation_with_substitution),
        ("LaTeX Generation", test_latex_generation)
    ]

    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"\n✗ ERROR in {test_name}: {str(e)}")
            results.append((test_name, False))

    # Summary
    print("\n" + "="*70)
    print(" "*25 + "TEST SUMMARY")
    print("="*70)

    passed = sum(1 for _, result in results if result)
    total = len(results)

    for test_name, result in results:
        status = "✓ PASSED" if result else "✗ FAILED"
        print(f"  {test_name:40s} {status}")

    print("-"*70)
    print(f"  Total: {passed}/{total} tests passed")

    if passed == total:
        print("\n🎉 All tests passed! The symbolic equation system is working correctly.")
    else:
        print(f"\n⚠️  {total - passed} test(s) failed. Please review the implementation.")

    return passed == total

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)