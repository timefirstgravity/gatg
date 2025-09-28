#!/usr/bin/env python
"""
Test suite for manifold and coordinate setup
Extracted from test_functions.py for better organization
"""

#!/usr/bin/env python
"""
Comprehensive pytest test suite for flux law mathematical functions
Tests the actual code from verify_flux_law_computational.py
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

class TestManifoldSetup:
    """Test manifold and coordinate setup"""

    def test_manifold_creation(self):
        """Test 4D Lorentzian manifold creation"""
        M, X, t, r, th, ph = setup_manifold_and_coordinates()

        # Test manifold dimension
        assert M.dim() == 4

        # Test Lorentzian structure is actually set
        assert hasattr(M, '_structure')
        assert str(type(M._structure).__name__) == 'LorentzianStructure'

        # Test chart belongs to manifold
        assert X.manifold() == M

        # Test coordinate system has 4 coordinates
        assert len(X[:]) == 4

        # Test that coordinates correspond to the chart coordinates
        assert t in X[:]
        assert r in X[:]
        assert th in X[:]
        assert ph in X[:]

    def test_coordinate_names(self):
        """Test coordinate variable names and ranges"""
        M, X, t, r, th, ph = setup_manifold_and_coordinates()

        # Test coordinate names
        coords = X[:]
        assert str(coords[0]) == 't'
        assert str(coords[1]) == 'r'
        assert str(coords[2]) == 'th'
        assert str(coords[3]) == 'ph'

        # Test that returned variables correspond to chart coordinates
        assert t == coords[0]
        assert r == coords[1]
        assert th == coords[2]
        assert ph == coords[3]

        # Test coordinate ranges match the spherical coordinate specification
        # Expected: t r:(0,+oo) th:(0,pi):theta ph:(0,2*pi):phi
        coord_range_str = str(X.coord_range())

        # Should contain the expected ranges
        assert 't: (-oo, +oo)' in coord_range_str  # t spans all reals
        assert 'r: (0, +oo)' in coord_range_str    # r > 0
        assert 'th: (0, pi)' in coord_range_str    # 0 < θ < π
        assert 'ph: (0, 2*pi)' in coord_range_str  # 0 < φ < 2π

        # Test individual coordinate range access
        assert str(X.coord_range(t)) == 't: (-oo, +oo)'
        assert str(X.coord_range(r)) == 'r: (0, +oo)'
        assert str(X.coord_range(th)) == 'th: (0, pi)'
        assert str(X.coord_range(ph)) == 'ph: (0, 2*pi)'

if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
