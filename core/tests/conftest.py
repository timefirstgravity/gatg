"""
Shared pytest fixtures and configuration for core tests.
"""

import pytest
import sys
import os

# Add sage to path and import functions
sys.path.append('/opt/homebrew/bin')  # Common sage location
os.environ['SAGE_ROOT'] = '/opt/homebrew'

try:
    from sage.all import *
    from core import *
except ImportError as e:
    pytest.skip(f"Sage not available: {e}", allow_module_level=True)

@pytest.fixture
def manifold_setup():
    """Fixture providing basic manifold setup used by multiple test classes."""
    M, X, t, r, th, ph = setup_manifold_and_coordinates()
    Phi, A, N = define_temporal_potential(M, t, r)
    return {
        'manifold': M,
        'coordinates': (t, r, th, ph),
        'temporal_potential': Phi,
        'metric_component': A,
        'lapse': N
    }

@pytest.fixture
def three_metric_setup():
    """Fixture providing 3-metric components for testing."""
    M, X, t, r, th, ph = setup_manifold_and_coordinates()
    Phi, A, N = define_temporal_potential(M, t, r)
    gamma_rr, gamma_thth, gamma_phph = compute_3metric_components(A, r, th)
    return {
        'A': A,
        'r': r,
        'th': th,
        'gamma_rr': gamma_rr,
        'gamma_thth': gamma_thth,
        'gamma_phph': gamma_phph
    }