from setuptools import setup, find_packages

setup(
    name="fit_abundance",
    version="0.1.0",
    description="Pipeline for oxygen abundance gradient fitting in galaxies",
    author="Andre Felipe de Siqueira Cardoso",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "statsmodels",
        "piecewise-regression",
    ],
    python_requires=">=3.8",
)
