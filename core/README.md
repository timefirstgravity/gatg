# GATG Core Module

The core module provides the fundamental mathematical framework for **Gravity as Temporal Geometry (GATG)** computations. This module contains the essential building blocks for symbolic computation and mathematical verification that Standard GR and Lapse-First GR formulations are equivalent.

## Core Modules

### Main Components

- **`__init__.py`** - Main module interface exposing all core functionality
- **`symbolic_equations.py`** - Central equation framework connecting string documentation with symbolic mathematics
- **`metrics.py`** - Metric tensor construction and temporal potential definitions
- **`tensors.py`** - Einstein tensor, Ricci tensor, and curvature computations
- **`manifolds.py`** - 4D Lorentzian manifold and coordinate system setup
- **`constraints.py`** - ADM Hamiltonian and momentum constraint computations
- **`adm.py`** - ADM 3+1 decomposition functions (extrinsic curvature, etc.)
- **`differential_operators.py`** - D'Alembertian, Laplacian, gradient, and divergence operators

### Equation Categories (`equations/` directory)

The symbolic equation framework is organized into specialized modules:

- **`schwarzschild.py`** - Schwarzschild solution equations and derivations
- **`flux_law.py`** - Fundamental flux law connecting temporal and spatial geometry
- **`adm.py`** - ADM formalism equations and constraints
- **`kerr_equations.py`** - Kerr black hole solution
- **`cosmological.py`** - Cosmological and FLRW equations
- **`christoffel.py`** - Christoffel symbols and connection computations
- **`conservation.py`** - Conservation laws and continuity equations
- **`coordinate_transformations.py`** - Coordinate transformation equations
- **`reissner_nordstrom.py`** - Reissner-Nordström charged black hole
- **`vaidya.py`** - Vaidya spacetime (null dust)
- **`kerr_newman.py`** - Kerr-Newman (rotating charged) solution
- **`dimensional_analysis.py`** - Dimensional analysis and unit verification
- **`cmb_observables.py`** - Cosmic microwave background observables

### Utility and Testing

- **`example_symbolic_usage.py`** - Comprehensive examples demonstrating framework usage
- **`tests/`** - Test suite for mathematical verification
- **Legacy files preserved for reference (marked `_old`)**

## Key Features

### 1. Rigorous Mathematical Computation

- All functions perform actual symbolic computation using SageMath
- No fake data or placeholder values - only real mathematical results
- Every computation contributes to GR equivalence verification

### 2. Modular Architecture

- Clean separation between different physics domains
- Reusable components across different spacetime solutions
- Consistent API across all modules

### 3. Symbolic-String Integration

- Connects documented equations with their symbolic representations
- Enables mathematical verification of documented formulas
- Automatic LaTeX generation from symbolic expressions

## Core API Functions

### Manifold Setup

```python
from core import setup_manifold_and_coordinates, setup_static_manifold_and_coordinates

# Standard 4D Lorentzian manifold with spherical coordinates
M, X, t, r, th, ph = setup_manifold_and_coordinates()

# Static version for time-independent problems
M, X, t, r, th, ph = setup_static_manifold_and_coordinates()
```

### Metric Construction

```python
from core import construct_metric, define_temporal_potential, compute_3metric_components

# Construct spherical metric tensor
g = construct_metric(M, A, r, th)

# Define temporal potential Φ and related quantities
Phi, A, N = define_temporal_potential(M, t, r)

# Compute 3-metric components
gamma_rr, gamma_thth, gamma_phph = compute_3metric_components(A, r, th)
```

### Tensor Computations

```python
from core import compute_einstein_tensor, compute_ricci_tensor, verify_vacuum_solution

# Compute Einstein tensor G_μν = R_μν - (1/2) R g_μν
G, R, Rs = compute_einstein_tensor(g)

# Extract components for analysis
components = extract_einstein_components(G, simplify=True)

# Verify vacuum solution (G_μν = 0)
is_vacuum, components = verify_vacuum_solution(G)
```

### ADM Formalism

```python
from core import compute_extrinsic_curvature, compute_hamiltonian_constraint

# Compute extrinsic curvature tensor
K_rr, K_thth, K_phph = compute_extrinsic_curvature(gamma_rr, gamma_thth, gamma_phph, N, t)

# Compute ADM constraints
hamiltonian_constraint = compute_hamiltonian_constraint(...)
momentum_constraint = compute_general_lambda_momentum_constraint(M, t, r)
```

### Differential Operators

```python
from core import dalembertian, spatial_laplacian, gradient_3d, divergence_3d

# D'Alembertian wave operator
wave_result = dalembertian(field, coordinates=[t,x,y,z])

# 3D spatial Laplacian
laplacian_result = spatial_laplacian(field, coordinates=[x,y,z])

# 3D gradient and divergence
grad_field = gradient_3d(field, coordinates=[x,y,z])
div_field = divergence_3d(vector_field, coordinates=[x,y,z])
```

### Symbolic Equation Framework

```python
from core import equations, get_string, get_symbolic

# Access equation by category and name
schwarzschild_eq = equations('schwarzschild', 'field_equation')

# Get string representation
eq_string = get_string('schwarzschild', 'field_equation')

# Get symbolic expression for computation
eq_symbolic = get_symbolic('schwarzschild', 'field_equation')
```

## Usage Requirements

### SageMath Environment

**CRITICAL:** All scripts must be run with `sage -python`, not plain `python`:

```bash
sage -python your_script.py
```

### Import Pattern

```python
#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

from sage.all import *
import core
from core import setup_manifold_and_coordinates, compute_einstein_tensor
```

## Example: Complete Schwarzschild Verification

```python
#!/usr/bin/env sage -python
from sage.all import *
from core import *

# Setup manifold and coordinates
M, X, t, r, th, ph = setup_manifold_and_coordinates()

# Define temporal potential (lapse-first approach)
Phi, A, N = define_temporal_potential(M, t, r, static=True)

# Apply Schwarzschild solution A(r) = 1 - 2M/r
M_mass = var('M', domain='positive')
A = 1 - 2*M_mass/r

# Construct metric with this A(r)
g = construct_metric(M, A, r, th)

# Compute Einstein tensor
G, R, Rs = compute_einstein_tensor(g)

# Verify vacuum solution G_μν = 0
is_vacuum, components = verify_vacuum_solution(G)
print(f"Vacuum solution verified: {is_vacuum}")

# This demonstrates that the complex Einstein equations
# reduce to the simple ODE: r*A'(r) + A(r) - 1 = 0
```
