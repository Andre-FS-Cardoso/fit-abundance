import numpy as np

def distances(ra,
              ra0,
              dec,
              dec0,
              pa,
              ba,
              d,
              re
              ):
    """
    Compute deprojected galactocentric distances of HII regions.

    The distances are calculated following the geometric deprojection
    described in Boczko (1998), Scarano (2008), Scarano et al. (2008),
    and Giovanelli et al. (1994). The resulting radius is normalized by
    the effective radius of the galaxy.

    Parameters
    ----------
    ra : array-like
        Right ascension of the HII regions (degrees).
    ra0 : float
        Right ascension of the galaxy center (degrees).
    dec : array-like
        Declination of the HII regions (degrees).
    dec0 : float
        Declination of the galaxy center (degrees).
    pa : float
        Position angle of the galaxy (degrees).
    ba : float
        Ratio between the minor and major axes (b/a).
    d : float
        Distance to the galaxy (Mpc).
    re : float
        Effective radius of the galaxy (kpc).

    Returns
    -------
    r : numpy.ndarray
        Deprojected galactocentric distances normalized by the
        effective radius (r / r_e).
    """
    
    # --- Convert angles to radians
    ra = np.asarray(ra) * np.pi / 180
    dec = np.asarray(dec) * np.pi / 180
    pa = pa*np.pi/180
    dec0 = dec0*np.pi/180
    ra0 = ra0*np.pi/180
    
    # --- Inclination correction
    cos_i = np.sqrt((ba**2-0.13**2)/(1-0.13**2))

    # --- Projected coordinates
    r1 = -(ra-ra0)*np.sin(pa)*np.cos(dec) + (dec-dec0)*np.cos(pa)
    r2 = (-(ra-ra0)*np.cos(pa)*np.cos(dec) - (dec-dec0)*np.sin(pa))/cos_i
    
    # --- Convert distance to kpc
    d_kpc = d * 1e3 # Mpc p/ kpc

    # --- Deprojected radius normalized by effective radius
    r = np.sqrt(r1**2 + r2**2) * d_kpc / re
    
    return r
