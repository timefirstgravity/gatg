#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
GATG Core: Differential Operators

Fundamental differential operators used throughout GATG computations:
- D'Alembertian (wave operator)
- Spatial Laplacian
- Covariant derivatives
- Divergence operators

These operators are used across multiple modules including linearized gravity,
cosmology, and exact solutions.
"""

from sage.all import *

def dalembertian(field, coordinates=None, signature='lorentzian'):
    """
    D'Alembertian wave operator: □ = η^μν ∂_μ ∂_ν

    For Minkowski spacetime: □ = -∂²_t + ∇²
    For general coordinates, uses the appropriate metric signature.

    Args:
        field: Symbolic expression to apply operator to
        coordinates: List of coordinate variables [t, x, y, z] (optional)
        signature: 'lorentzian' (-,+,+,+) or 'euclidean' (+,+,+,+)

    Returns:
        Result of applying □ to the field
    """
    if coordinates is None:
        t, x, y, z = var('t x y z')
    else:
        t, x, y, z = coordinates[:4]  # Take first 4 coordinates

    if signature == 'lorentzian':
        # Minkowski signature (-,+,+,+)
        return -diff(field, t, t) + diff(field, x, x) + diff(field, y, y) + diff(field, z, z)
    elif signature == 'euclidean':
        # Euclidean signature (+,+,+,+)
        return diff(field, t, t) + diff(field, x, x) + diff(field, y, y) + diff(field, z, z)
    else:
        raise ValueError(f"Unknown signature: {signature}. Use 'lorentzian' or 'euclidean'")

def spatial_laplacian(field, coordinates=None):
    """
    Spatial Laplacian: ∇² = ∂²/∂x² + ∂²/∂y² + ∂²/∂z²

    Args:
        field: Symbolic expression to apply operator to
        coordinates: List of spatial coordinate variables [x, y, z] (optional)

    Returns:
        Result of applying ∇² to the field
    """
    if coordinates is None:
        x, y, z = var('x y z')
    else:
        # Handle both 3D [x,y,z] and 4D [t,x,y,z] coordinate lists
        if len(coordinates) == 3:
            x, y, z = coordinates
        else:
            x, y, z = coordinates[1:4]  # Skip time coordinate

    return diff(field, x, x) + diff(field, y, y) + diff(field, z, z)

def gradient_3d(field, coordinates=None):
    """
    3D gradient: ∇f = (∂f/∂x, ∂f/∂y, ∂f/∂z)

    Args:
        field: Scalar field to take gradient of
        coordinates: List of spatial coordinate variables [x, y, z] (optional)

    Returns:
        List [∂f/∂x, ∂f/∂y, ∂f/∂z]
    """
    if coordinates is None:
        x, y, z = var('x y z')
    else:
        if len(coordinates) == 3:
            x, y, z = coordinates
        else:
            x, y, z = coordinates[1:4]  # Skip time coordinate

    return [diff(field, x), diff(field, y), diff(field, z)]

def divergence_3d(vector_field, coordinates=None):
    """
    3D divergence: ∇·V = ∂V^x/∂x + ∂V^y/∂y + ∂V^z/∂z

    Args:
        vector_field: List [V_x, V_y, V_z] of vector components
        coordinates: List of spatial coordinate variables [x, y, z] (optional)

    Returns:
        Scalar divergence ∇·V
    """
    if coordinates is None:
        x, y, z = var('x y z')
    else:
        if len(coordinates) == 3:
            x, y, z = coordinates
        else:
            x, y, z = coordinates[1:4]  # Skip time coordinate

    if len(vector_field) != 3:
        raise ValueError("Vector field must have 3 components [V_x, V_y, V_z]")

    V_x, V_y, V_z = vector_field
    return diff(V_x, x) + diff(V_y, y) + diff(V_z, z)

def time_derivative(field, time_coord=None, order=1):
    """
    Time derivative: ∂^n f/∂t^n

    Args:
        field: Field to differentiate
        time_coord: Time coordinate variable (optional, defaults to 't')
        order: Order of derivative (default 1)

    Returns:
        Time derivative of specified order
    """
    if time_coord is None:
        t = var('t')
    else:
        t = time_coord

    result = field
    for _ in range(order):
        result = diff(result, t)
    return result