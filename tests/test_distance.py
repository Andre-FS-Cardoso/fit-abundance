import pytest
import numpy as np
from distance import distances

def test_center_is_zero():
    """Checks whether the galaxy center results in r = 0"""
    
    # Mock data for a galaxy
    ra0, dec0 = 10.0, -10.0
    pa, ba = 45.0, 0.5
    dist_mpc, re_kpc = 100.0, 5.0
    
    # HII region exactly at the galaxy center
    r = distances([ra0], ra0, [dec0], dec0, pa, ba, dist_mpc, re_kpc)
    
    assert np.isclose(r[0], 0.0)

def test_symmetry_face_on():
    """Checks whether symmetric regions in a face-on galaxy (ba = 1) have the same r"""
     
    ra0, dec0 = 0.0, 0.0
    pa, ba = 0.0, 1.0
    dist, re = 50.0, 2.0
    
    # Two symmetric regions (one east and one west of the center)
    ra_pts = [0.1, -0.1]
    dec_pts = [0.0, 0.0]
    
    r = distances(ra_pts, ra0, dec_pts, dec0, pa, ba, dist, re)
    
    # The result should be the same for both (numerical tolerance)
    assert np.isclose(r[0], r[1], atol=1e-5)

def test_known_value_ngc0309():
    """Checks whether the calculation for NGC 0309 produces a reasonable value (sanity test)"""
    
    # Using the values from the data_NGC0309.csv file
    ra0, dec0 = 14.17776, -9.91387
    pa, ba = 107.84957, 0.84
    dist, re = 79.55, 13.47
    
    # Test a region offset only in RA
    ra_hii = [14.18] 
    dec_hii = [-9.91387]
    
    r = distances(ra_hii, ra0, dec_hii, dec0, pa, ba, dist, re)
    
    # r must be positive and non-zero
    assert r[0] > 0
    assert isinstance(r, np.ndarray)
