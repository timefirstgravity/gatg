# -*- coding: utf-8 -*-
"""
GATG Equations Package
Modular organization of symbolic equations by category
"""

from .schwarzschild import schwarzschildEquations
from .flux_law import FluxLawEquations
from .adm import AdmEquations
from .kerr_equations import KerrEquations
from .cosmological import CosmologicalEquations
from .christoffel import ChristoffelEquations
from .conservation import ConservationEquations
from .coordinate_transformations import CoordinateTransformationEquations
from .reissner_nordstrom import ReissnerNordstromEquations
from .vaidya import VaidyaEquations
from .kerr_newman import KerrNewmanEquations
from .dimensional_analysis import DimensionalAnalysisEquations
from .cmb_observables import CmbObservablesEquations

__all__ = [
    'schwarzschildEquations',
    'FluxLawEquations',
    'AdmEquations',
    'KerrEquations',
    'CosmologicalEquations',
    'ChristoffelEquations',
    'ConservationEquations',
    'CoordinateTransformationEquations',
    'ReissnerNordstromEquations',
    'VaidyaEquations',
    'KerrNewmanEquations',
    'DimensionalAnalysisEquations',
    'CmbObservablesEquations'
]
