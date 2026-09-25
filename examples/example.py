from fit_abundance import fit_final
import pandas as pd

###############################################################################

galaxy_data = pd.read_csv("data_NGC0309.csv")

name = galaxy_data["galaxy"].iloc[0]
ra0 = galaxy_data["ra0"].iloc[0]
dec0 = galaxy_data["dec0"].iloc[0]
pa = galaxy_data["pa"].iloc[0]
ba = galaxy_data["ba"].iloc[0]
re = galaxy_data["re"].iloc[0]
d = galaxy_data["dist"].iloc[0]

###############################################################################

flux = pd.read_csv("HII.NGC0309.flux_elines.csv")

ra = flux["RA"]
dec = flux["DEC"]
EWHa = flux["EWHa6562"]
HIIREGID = flux["HIIREGID"]
Hb4861 = flux["fluxHb4861"]
eHb4861 = flux["e_fluxHb4861"]
OII3727 = flux["fluxOII3727"]
eOII3727 = flux["e_fluxOII3727"]
OIII4958 = flux["fluxOIII4958"]
eOIII4958 = flux["e_fluxOIII4958"]
OIII5006 = flux["fluxOIII5006"]
eOIII5006 = flux["e_fluxOIII5006"]
Ha6562 = flux["fluxHa6562"]
eHa6562 = flux["e_fluxHa6562"]
NII6548 = flux["fluxNII6548"]
eNII6548  = flux["e_fluxNII6548"]
NII6583 = flux["fluxNII6583"]
eNII6583 = flux["e_fluxNII6583"]
SII6716 = flux["fluxSII6716"]
eSII6716 = flux["e_fluxSII6716"]
SII6730 = flux["fluxSII6730"]
eSII6730 = flux["e_fluxSII6730"]

###############################################################################

calibrator = 1
criterion = "KA03"

###############################################################################

results = fit_final(
        name, HIIREGID, ra, ra0, dec, dec0, pa, ba, d, re, EWHa,
        Hb4861, eHb4861, Ha6562, eHa6562, OII3727, eOII3727,
        OIII4958, eOIII4958, OIII5006, eOIII5006, NII6548, eNII6548,
        NII6583, eNII6583, SII6716, eSII6716, SII6730, eSII6730,
        calibrator, criterion,
        n_boot=200,
        n_break=2,
        show_graph_bpt=True,
        save_abundance=True,
        save_oh_criteria=True,
        save_BPT=True,
        save_model_selection=True,
        save_graph=True,
        show_graph=True
        )

print(results)
