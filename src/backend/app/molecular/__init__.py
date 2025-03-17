"""
Molecular Data Management and Processing Module

This module provides a unified interface to the molecular processing capabilities
of the Molecular Data Management and CRO Integration Platform. It includes functionality
for molecular structure handling, property calculation, validation, similarity searching,
and substructure matching.

Key components:
- Molecule conversion and representation
- Structure validation and property checking
- Batch processing of molecular data
- Molecular property calculation
- Similarity searching and molecule clustering
- Substructure searching and highlighting
"""

import logging

# Import key functions and classes from submodules
from .molecule_converter import (
    smiles_to_mol, 
    smiles_to_mol_with_exception,
    mol_to_smiles, 
    canonicalize_smiles,
    mol_to_image, 
    smiles_to_image, 
    smiles_to_base64_image,
    mol_to_fingerprint,
    FINGERPRINT_TYPES
)

from .validator import (
    validate_smiles_structure,
    validate_smiles_structure_with_exception,
    validate_molecule,
    validate_molecules_batch,
    PROPERTY_RANGES,
    check_lipinski_rule_of_five,
    check_veber_rules
)

from .processor import (
    process_molecule,
    process_molecules_batch,
    process_molecules_stream,
    MoleculeProcessor,
    filter_molecules,
    sort_molecules,
    deduplicate_molecules
)

from .property_calculator import (
    calculate_properties,
    calculate_properties_from_smiles,
    batch_calculate_properties,
    PropertyCalculator,
    PROPERTY_CALCULATORS
)

from .similarity_searcher import (
    calculate_similarity,
    calculate_similarity_from_smiles,
    find_similar_molecules,
    SimilaritySearcher
)

from .substructure_searcher import (
    search_substructure,
    search_substructure_smiles,
    filter_molecules_with_substructure,
    SubstructureSearcher
)

# Set up module logger
logger = logging.getLogger(__name__)

# Module version
__version__ = '1.0.0'