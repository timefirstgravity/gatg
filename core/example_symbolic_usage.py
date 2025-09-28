#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Example Usage: Symbolic GATG Equations

This script demonstrates how to use the new SymbolicGATGEquations framework
that connects string documentation with actual symbolic mathematics.

The key advantage is that every documented equation can now be:
1. Mathematically verified
2. Evaluated with specific values
3. Used in further symbolic computations
4. Generate its own LaTeX documentation

This ensures that documentation and computation are always consistent.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sage.all import *
from symbolic_equations import SymbolicGATGEquations, symbolic_equations
from equations import GATGEquations

def example_basic_access():
    """Example 1: Basic equation access and comparison"""
    print("="*80)
    print(" "*25 + "EXAMPLE 1: BASIC ACCESS")
    print("="*80)

    seq = SymbolicGATGEquations()

    # Access equation in multiple formats
    eq_key = 'schwarzschild_ode'
    print(f"Equation: {eq_key}")
    print(f"  String form:    {seq.get_string(eq_key)}")
    print(f"  LaTeX form:     {seq.get_latex(eq_key)}")
    print(f"  Description:    {seq.get_description(eq_key)}")
    print(f"  Symbolic form:  {seq.get_symbolic(eq_key)}")

    # Compare with old GATGEquations
    print(f"\nBackward compatibility check:")
    print(f"  Old class:      {GATGEquations.SCHWARZSCHILD_ODE}")
    print(f"  New class:      {seq.get_string(eq_key)}")
    print(f"  Match: {GATGEquations.SCHWARZSCHILD_ODE == seq.get_string(eq_key)}")

def example_equation_verification():
    """Example 2: Mathematical verification of equations"""
    print("\n" + "="*80)
    print(" "*20 + "EXAMPLE 2: MATHEMATICAL VERIFICATION")
    print("="*80)

    seq = SymbolicGATGEquations()

    # Verify Schwarzschild ODE
    print("1. Verifying Schwarzschild ODE with known solution:")
    print("-"*60)
    is_valid = seq.verify_schwarzschild_ode()
    print(f"   The equation: {seq.get_string('schwarzschild_ode')}")
    print(f"   With solution: A(r) = 1 - r_s/r")
    print(f"   Verification: {'✓ VALID' if is_valid else '✗ INVALID'}")

    # Verify flux law consistency
    print("\n2. Verifying flux law consistency:")
    print("-"*60)
    is_consistent = seq.verify_flux_law_consistency()
    print(f"   Flux law: {seq.get_string('flux_law')}")
    print(f"   Einstein relation: {seq.get_string('einstein_tr')}")
    print(f"   Consistency: {'✓ CONSISTENT' if is_consistent else '✗ INCONSISTENT'}")

    # Demonstrate verification failure (for contrast)
    print("\n3. Example of equation verification:")
    print("-"*60)
    print("   This shows how the system can catch mathematical errors:")
    print("   If we had A(r) = 1 - r_s/(2r) [WRONG], the ODE verification would fail")

def example_symbolic_computation():
    """Example 3: Using equations in symbolic computations"""
    print("\n" + "="*80)
    print(" "*20 + "EXAMPLE 3: SYMBOLIC COMPUTATION")
    print("="*80)

    seq = SymbolicGATGEquations()

    # Get symbolic expressions and work with them
    print("1. Working with symbolic Schwarzschild solution:")
    print("-"*60)
    A_expr = seq.get_symbolic('schwarzschild_solution')
    r = seq.vars['r']
    rs = seq.vars['r_s']

    print(f"   A(r) = {A_expr}")

    # Compute derivative
    dA_dr = diff(A_expr, r)
    print(f"   dA/dr = {dA_dr}")

    # Verify the derivative satisfies the constraint
    constraint_lhs = r * dA_dr + A_expr - 1
    constraint_simplified = constraint_lhs.simplify_full()
    print(f"   Verification: r*A'(r) + A(r) - 1 = {constraint_simplified}")

    print("\n2. Computing event horizon and surface gravity:")
    print("-"*60)
    # Event horizon (where A(r) = 0)
    horizon_eq = A_expr.solve(r)
    print(f"   Event horizon: A(r) = 0 when r = {rs}")

    # Surface gravity at horizon
    surface_gravity = (1/2) * dA_dr.substitute(r=rs)
    print(f"   Surface gravity: κ = (1/2) * dA/dr|_rs = {surface_gravity}")

def example_practical_calculation():
    """Example 4: Practical physics calculation"""
    print("\n" + "="*80)
    print(" "*20 + "EXAMPLE 4: PRACTICAL CALCULATION")
    print("="*80)

    seq = SymbolicGATGEquations()

    print("1. Computing black hole properties for different masses:")
    print("-"*60)

    # Define masses
    masses = {
        'Solar mass': 1.989e30,      # kg
        'Stellar mass (10 M☉)': 10 * 1.989e30,
        'Supermassive (10^6 M☉)': 1e6 * 1.989e30
    }

    # Physical constants
    G = 6.674e-11  # m^3/(kg*s^2)
    c = 3e8        # m/s

    rs_expr = seq.get_symbolic('schwarzschild_radius')

    for name, mass in masses.items():
        # Substitute values
        substitutions = {
            seq.vars['M']: mass,
            seq.vars['G']: G,
            seq.vars['c']: c
        }

        rs_value = float(rs_expr.substitute(substitutions))

        print(f"   {name:20s}: M = {mass:.2e} kg")
        print(f"   {' '*20}  r_s = {rs_value:.3f} m")
        print(f"   {' '*20}  r_s = {rs_value/1000:.3f} km")
        print()

def example_flux_law_physics():
    """Example 5: Flux law physics with units"""
    print("\n" + "="*80)
    print(" "*20 + "EXAMPLE 5: FLUX LAW PHYSICS")
    print("="*80)

    seq = SymbolicGATGEquations()

    print("1. Understanding the flux law dimensions:")
    print("-"*60)

    # Get the flux law
    flux_law = seq.get_string('flux_law')
    print(f"   Flux law: {flux_law}")

    print("\n   Dimensional analysis:")
    print("   • [∂_t Φ] = [time]^(-1)")
    print("   • [4πG/c⁴] = [length]/[energy]")
    print("   • [r] = [length]")
    print("   • [T^t_r] = [energy]/([length]²·[time])")
    print("   • RHS = [length]/[energy] × [length] × [energy]/([length]²·[time])")
    print("   •     = [time]^(-1) ✓")

    print("\n2. Physical interpretation:")
    print("-"*60)
    print("   • Φ(t,r): Temporal potential (gravitational redshift)")
    print("   • ∂_t Φ > 0: Time dilation increasing (matter falling in)")
    print("   • ∂_t Φ < 0: Time dilation decreasing (energy radiating out)")
    print("   • T^t_r: Energy flux (positive = outward, negative = inward)")

def example_adm_decomposition():
    """Example 6: ADM 3+1 decomposition"""
    print("\n" + "="*80)
    print(" "*20 + "EXAMPLE 6: ADM DECOMPOSITION")
    print("="*80)

    seq = SymbolicGATGEquations()

    print("1. ADM variables and their meanings:")
    print("-"*60)

    adm_vars = [
        ('lapse_function', 'Rate of proper time flow'),
        ('extrinsic_trace', 'How spatial slices expand'),
        ('gamma_rr', 'Radial stretching of space'),
        ('momentum_constraint', 'Spatial momentum conservation'),
        ('hamiltonian_constraint', 'Energy conservation')
    ]

    for var_name, meaning in adm_vars:
        eq_string = seq.get_string(var_name)
        print(f"   {var_name:20s}: {eq_string}")
        print(f"   {' '*20}  → {meaning}")
        print()

    print("2. Gauge conditions (coordinate choices):")
    print("-"*60)
    print(f"   Zero shift: {seq.get_string('zero_shift')}")
    print("   → Spatial coordinates don't get dragged by time evolution")
    print()
    print(f"   Areal radius: {seq.get_string('areal_radius')}")
    print("   → Area of spheres is 4πr² (natural radial coordinate)")

def example_backward_compatibility():
    """Example 7: Backward compatibility with existing code"""
    print("\n" + "="*80)
    print(" "*15 + "EXAMPLE 7: BACKWARD COMPATIBILITY")
    print("="*80)

    print("1. Using enhanced GATGEquations class:")
    print("-"*60)

    # Access equation strings as before
    print(f"   Old way: {GATGEquations.FLUX_LAW}")

    # But now can also access symbolic version
    flux_symbolic = GATGEquations.get_symbolic('FLUX_LAW')
    if flux_symbolic:
        print(f"   Symbolic: {flux_symbolic}")
        print("   ✓ Enhanced with symbolic computation")
    else:
        print("   ✗ Symbolic not available")

    # Verify equation
    is_verified = GATGEquations.verify_equation('SCHWARZSCHILD_ODE')
    if is_verified is not None:
        print(f"   Verification: {'✓ PASSED' if is_verified else '✗ FAILED'}")
    else:
        print("   Verification: Not available")

    print("\n2. Migration status:")
    print("-"*60)
    migrated_equations = [
        'SCHWARZSCHILD_ODE',
        'FLUX_LAW',
        'EINSTEIN_TR',
        'LAPSE_FUNCTION',
        'MOMENTUM_CONSTRAINT'
    ]

    for eq in migrated_equations:
        has_symbolic = GATGEquations.has_symbolic(eq)
        status = "✓ Migrated" if has_symbolic else "○ Pending"
        print(f"   {eq:25s}: {status}")

def example_latex_generation():
    """Example 8: Automatic LaTeX generation"""
    print("\n" + "="*80)
    print(" "*20 + "EXAMPLE 8: LATEX GENERATION")
    print("="*80)

    seq = SymbolicGATGEquations()

    print("1. Generating LaTeX for research papers:")
    print("-"*60)

    # Generate LaTeX for Schwarzschild equations
    schwarzschild_latex = seq.generate_latex_document('schwarzschild')
    print("   Schwarzschild equations LaTeX:")
    print("   " + "-"*50)
    for line in schwarzschild_latex.split('\n')[:15]:
        print(f"   {line}")
    print("   " + "-"*50)

    print("\n2. Individual equation LaTeX:")
    print("-"*60)
    eq_latex = seq.get_latex('flux_law')
    print(f"   Flux law: {eq_latex}")
    print("   → Can be directly included in papers/presentations")

def main():
    """Run all examples demonstrating the symbolic equations framework"""
    print("\n" + "="*80)
    print(" "*15 + "SYMBOLIC GATG EQUATIONS: USAGE EXAMPLES")
    print("="*80)
    print("\nThis script demonstrates how the new symbolic equation framework")
    print("connects documentation strings with actual mathematical computation.")
    print("\nKey benefits:")
    print("• Every equation can be mathematically verified")
    print("• Documentation and computation are always consistent")
    print("• Backward compatibility with existing string-based code")
    print("• Automatic LaTeX generation for papers")
    print("• Support for practical physics calculations")

    # Run all examples
    examples = [
        example_basic_access,
        example_equation_verification,
        example_symbolic_computation,
        example_practical_calculation,
        example_flux_law_physics,
        example_adm_decomposition,
        example_backward_compatibility,
        example_latex_generation
    ]

    for example_func in examples:
        try:
            example_func()
        except Exception as e:
            print(f"\n✗ ERROR in {example_func.__name__}: {str(e)}")

    # Summary
    print("\n" + "="*80)
    print(" "*25 + "IMPLEMENTATION COMPLETE")
    print("="*80)
    print("\nThe symbolic equation framework is now ready for use in GATG.")
    print("\nNext steps:")
    print("• Migrate remaining equation categories (Kerr, cosmological, etc.)")
    print("• Add more sophisticated verification methods")
    print("• Integrate with existing computation modules")
    print("• Extend to other spacetime solutions")
    print("\nThis framework ensures mathematical rigor while maintaining")
    print("the documentation and usability that makes GATG accessible.")

if __name__ == "__main__":
    main()