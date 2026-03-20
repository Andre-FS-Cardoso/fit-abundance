import pytest
import numpy as np
from models import fit_models

def test_select_linear_model():
    """Generates linear data with noise and checks whether AIC selects the simple model"""
    
    np.random.seed(42)
    x = np.linspace(0.1, 2.5, 40)
    
    # Line: y = -0.05x + 8.6 + Gaussian noise
    noise = np.random.normal(0, 0.01, size=len(x))
    y = -0.05 * x + 8.6 + noise
    ey = np.full_like(x, 0.01) 
    
    results = fit_models(x, y, ey, "TestLinear", "ST06", 1, False)
    
    # The linear model should clearly be selected
    assert results['best_case'] == 1

def test_select_broken_model():
    """Generates data with an abrupt break at radius r = 1.0"""
    
    x = np.linspace(0.1, 2.5, 40)
    
    # Strong slope until r=1, then flat
    y = np.where(x < 1.0, 8.7 - 0.2*x, 8.5)
    ey = np.full_like(x, 0.005)
    
    results = fit_models(x, y, ey, "TestBroken", "ST06", 1, False)
    
    # The algorithm should detect that the linear model (1) is not the best
    assert results['best_case'] >= 2

def test_minimum_points_requirement():
    """Checks whether the code ignores galaxies with fewer than 10 data points"""
    
    x = np.linspace(0.1, 0.5, 5)
    y = -0.05 * x + 8.5
    ey = np.full_like(x, 0.01)
    
    results = fit_models(x, y, ey, "SmallSample", "ST06", 1, False)
    
    # According to models.py: if len(x) < 10: return None
    assert results is None
