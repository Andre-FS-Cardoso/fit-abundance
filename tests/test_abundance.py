import pytest
from fit_abundance.abundance import abundance

def test_calibrator_values():
    """Checks whether the calibrators return physically plausible values"""
    
    # Test data: a hypothetical HII region with typical fluxes
    # Halpha/Hbeta ~ 2.86 (theoretical), relatively high NII/Halpha (metal-rich)
    name = "TestGal"
    r = [0.5]
    HIIREGID = ["HII-1"]
    EWHa = [20.0]

    Hb = [100.0]
    eHb = [1.0]

    Ha = [286.0]
    eHa = [2.0]

    OII3727 = [100.0]
    eOII3727 = [1.0]

    OIII4958 = [17.0]
    eOIII4958 = [1.0]

    OIII5006 = [50.0]
    eOIII5006 = [1.0]

    NII6548 = [33.0]
    eNII6548 = [1.0]

    NII6583 = [100.0]
    eNII6583 = [1.0]

    SII6716 = [30.0]
    eSII6716 = [1.0]

    SII6730 = [20.0]
    eSII6730 = [1.0]

    # Testing PP04 O3N2 (Calibrator 1)
    oh, eoh, ha_c, oiii_c, nii_c = abundance(
        name, r, HIIREGID, EWHa,
        Hb, eHb,
        Ha, eHa, 
        OII3727, eOII3727,
        OIII4958, eOIII4958,
        OIII5006, eOIII5006,
        NII6548, eNII6548,
        NII6583, eNII6583,
        SII6716, eSII6716, 
        SII6730, eSII6730,
        calibrator=1,
        save_abundance=False
        )
    
    # Oxygen abundance typically lies between 7.0 and 9.5
    assert 7.0 <= oh[0] <= 9.5
    assert eoh[0] > 0

def test_invalid_calibrator():
    """Checks whether the code raises an error for a non-existent calibrator"""
    
    with pytest.raises(ValueError, match="Invalid calibrator"):
        abundance(
            "Test",
            [0],
            ["1"],
            [1],
            [1], [1],
            [1], [1],
            [1], [1],
            [1], [1],
            [1], [1],
            [1], [1],
            [1], [1],
            [1], [1],
            [1], [1],
            calibrator=99,
            save_abundance=False
            )	
