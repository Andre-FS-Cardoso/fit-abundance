# fit_abundance

**fit_abundance** is a Python package designed to compute oxygen abundances and fit radial abundance gradients in galaxies using spectroscopic data from HII regions.

The software provides a fully automated pipeline for abundance-gradient analysis commonly used in extragalactic astronomy.

---

## Features

The pipeline performs the following steps:

1. Deprojection of HII region positions
2. Extinction correction of emission-line fluxes
3. Oxygen abundance calculations using strong-line calibrators
4. Selection of HII regions based on spectral criteria
5. Fitting of abundance gradients
6. Automatic model selection using the Akaike Information Criterion (AIC)

---

## Implemented abundance calibrators

The following strong-line calibrators are implemented:

- **PP04 (O3N2)** — Pettini & Pagel (2004)
- **PP04 (N2)** — Pettini & Pagel (2004)
- **M13 (O3N2)** — Marino et al. (2013)
- **M13 (N2)** — Marino et al. (2013)
- **D16** — Dopita et al. (2016)

---

## Implemented selection criteria for H II regions

- **ST06** : Stasińska et al. (2006)
- **KA03** : Kauffmann et al. (2003)
- **KE01** : Kewley et al. (2001)
- **KE6A** : Kewley et al. (2001) with EW(Hα) ≥ 6 Å
- **CF11** : Cid Fernandes et al. (2011)
- **EP20** : Espinosa-Ponce et al. (2020)

---

## Testing

To ensure scientific integrity and code reliability, we provide a suite of automated tests.

1- Install pytest:

```bash
pip install pytest
```

2- Run the tests:

```bash
pytest
```

---

## Installation

Clone the repository:

```bash
git clone https://github.com/Andre-FS-Cardoso/fit-abundance.git
cd fit-abundance
```

Install the required dependencies:

```bash
pip install -r requirements.txt
```

Install the package locally:

```bash
pip install -e .
```

---

## Dependencies

The package requires the following Python libraries:

- numpy
- pandas
- matplotlib
- statsmodels
- piecewise-regression

---

## Example

```python
from fit_abundance import fit_final

results = fit_final(...)
print(results)
```

```markdown
A complete working example is available in:

examples/example.py
```

---

## Output

The software returns a dictionary containing the fitted parameters of the oxygen abundance gradient.

Optionally, it can generate:

- tables of selected HII regions
- tables with model-selection statistics
- plots of the abundance gradient

---

## Authors

André Felipe de Siqueira Cardoso

Oscar Cavichia

---

## Citation

If you use this software in scientific work, please cite:

Cardoso, A. F. S. and Cavichia, O. (2026)

*fit_abundance: A Python tool for fitting oxygen abundance gradients in galaxies*

---

## License

This project is licensed under the MIT License.
