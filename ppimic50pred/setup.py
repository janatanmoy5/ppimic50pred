from setuptools import setup, find_packages

setup(
    name="ppimic50pred",
    version="0.1.0",
    author="Tanmoy Jana",
    description="Small Molecule PPI IC50 Predictor",
    packages=find_packages(),
    include_package_data=True,  # To include model file if specified in MANIFEST.in
    install_requires=[
        "pandas",
        "requests",
        "rdkit-pypi",
        "scikit-learn",
        "chembl-webresource-client"
    ],
    python_requires='>=3.7',
)
