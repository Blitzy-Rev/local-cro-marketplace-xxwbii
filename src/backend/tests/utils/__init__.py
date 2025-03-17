"""
Test utilities package for the Molecular Data Management and CRO Integration Platform.

This package provides common test utilities and fixtures for reuse across test modules.
"""

# List of valid SMILES strings for testing molecular functions
VALID_TEST_SMILES = [
    'CCO',              # Ethanol
    'c1ccccc1',         # Benzene
    'CC(=O)O',          # Acetic acid
    'C1CCCCC1',         # Cyclohexane
    'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'  # Caffeine
]

# List of invalid SMILES strings for testing error handling in molecular functions
INVALID_TEST_SMILES = [
    '',                 # Empty string
    'X',                # Invalid atom
    'C1CC',             # Unclosed ring
    'c1ccc'             # Incomplete aromatic ring
]

# List of common molecular properties for testing property calculation functions
TEST_MOLECULE_PROPERTIES = [
    'MolecularWeight',
    'LogP',
    'NumHDonors',
    'NumHAcceptors',
    'TPSA'
]