# docs/molecular-data-guidelines.md

# Introduction

This document provides comprehensive guidelines for handling molecular data in the Molecular Data Management and CRO Integration Platform. It covers SMILES notation standards, molecular property definitions, CSV import formats, and data validation rules. Following these guidelines ensures efficient data processing, accurate validation, and seamless integration with CRO services.

## Purpose

The purpose of this document is to:

- Establish standards for molecular data representation
- Define validation rules for molecular structures and properties
- Provide guidance for preparing CSV files for import
- Document troubleshooting procedures for common data issues

## Scope

These guidelines apply to all molecular data handled by the platform, including:

- Imported CSV files containing molecular structures and properties
- Molecular libraries created and managed within the platform
- Experimental data submitted to CROs
- Results received from CRO services

# SMILES Notation Guidelines

## SMILES Format

SMILES (Simplified Molecular Input Line Entry System) is a string notation used to represent molecular structures. The platform uses RDKit // rdkit 2023.03+ for SMILES parsing and validation.

## Validation Rules

Valid SMILES strings must:

- Contain valid atom symbols (e.g., C, N, O, S, P, etc.)
- Have balanced parentheses and brackets
- Follow correct bond notation (single, double, triple, aromatic)
- Represent a chemically valid structure

## Examples

| Molecule | Valid SMILES | Notes |
|----------|--------------|-------|
| Ethanol | CCO | Valid simple structure |
| Benzene | c1ccccc1 | Valid aromatic structure |
| Aspirin | CC(=O)Oc1ccccc1C(=O)O | Valid complex structure |
| C(C | Invalid | Unbalanced parentheses |
| C1CC | Invalid | Unclosed ring |

## Canonicalization

The platform automatically canonicalizes SMILES strings to ensure consistent representation. This process generates a unique SMILES string for each molecular structure, enabling accurate deduplication and comparison.

# Molecular Properties

## Supported Properties

The platform supports the following molecular properties:

| Property | Description | Valid Range | Units |
|----------|-------------|-------------|-------|
| molecular_weight | Molecular weight | 0-1000 | g/mol |
| logp | Octanol-water partition coefficient | -10-10 | - |
| h_bond_donors | Number of hydrogen bond donors | 0-10 | count |
| h_bond_acceptors | Number of hydrogen bond acceptors | 0-20 | count |
| rotatable_bonds | Number of rotatable bonds | 0-15 | count |
| polar_surface_area | Topological polar surface area | 0-200 | Å² |
| heavy_atom_count | Number of non-hydrogen atoms | 0-100 | count |
| ring_count | Number of rings | 0-10 | count |
| aromatic_rings | Number of aromatic rings | 0-7 | count |
| solubility | Aqueous solubility (log scale) | -10-10 | log(mol/L) |

## Property Calculation

The platform uses RDKit // rdkit 2023.03+ to calculate molecular properties. Properties can be:

1. Imported from CSV files
2. Calculated automatically for missing properties
3. Recalculated on demand

## Drug-likeness Rules

### Lipinski's Rule of Five

The platform evaluates molecules against Lipinski's Rule of Five for drug-likeness:

- Molecular weight ≤ 500 Da
- LogP ≤ 5
- Hydrogen bond donors ≤ 5
- Hydrogen bond acceptors ≤ 10

Molecules with no more than one violation are considered drug-like.

### Veber Rules

The platform also evaluates molecules against Veber's rules for oral bioavailability:

- Rotatable bonds ≤ 10
- Polar surface area ≤ 140 Å²

# CSV Import Format

## File Requirements

- Format: Comma-separated values (CSV)
- Encoding: UTF-8
- Maximum file size: 50 MB
- Maximum rows: No hard limit, but performance optimized for up to 10,000 molecules

## Required Columns

At minimum, CSV files must contain a column with SMILES strings. This column can have any header name but must be mapped to the SMILES property during import.

## Recommended Columns

The following additional columns are recommended:

- Molecule name or identifier
- Source or reference
- Custom properties relevant to your research

## Example CSV Format

```
Compound_ID,SMILES,MW,LogP,Activity
CPD001,CCO,46.07,-0.14,78.5
CPD002,CCCCO,74.12,0.88,45.2
CPD003,c1ccccc1,78.11,1.90,92.1
```

## Column Mapping

During import, you will map CSV columns to system properties:

1. SMILES column must be mapped (required)
2. Other columns can be mapped to system properties or custom properties
3. Unmapped columns will be ignored

## Batch Processing

Large CSV files are processed in batches of 1,000 rows to optimize performance and memory usage. Progress is reported during processing.

# Data Validation Rules

## Structure Validation

All molecular structures undergo validation during import:

1. SMILES syntax validation
2. Chemical validity check
3. Structure normalization

Invalid structures are flagged and reported but not imported.

## Property Validation

Molecular properties are validated against expected ranges:

1. Numerical properties must be within valid ranges (see Molecular Properties section)
2. Missing properties are calculated when possible
3. Properties outside valid ranges are flagged but still imported

## Duplicate Handling

Duplicate molecules are identified based on canonical SMILES:

1. Exact duplicates (same SMILES) within a single import are removed
2. Duplicates across different imports are maintained as separate entries
3. Duplicate detection can be enabled for cross-import comparison

## Validation Report

After import, a validation report is generated with:

1. Total molecules processed
2. Valid molecules imported
3. Invalid molecules rejected
4. Duplicate molecules detected
5. Property validation issues
6. Detailed error messages for troubleshooting

# Best Practices

## Data Preparation

1. **Standardize SMILES**: Use canonical SMILES when possible
2. **Validate Before Import**: Pre-validate molecular structures using RDKit // rdkit 2023.03+ or similar tools
3. **Include Identifiers**: Add unique identifiers for each molecule
4. **Document Properties**: Include metadata about property sources and calculation methods
5. **Batch Appropriately**: Split very large datasets into multiple files of 5,000-10,000 molecules

## Data Organization

1. **Create Focused Libraries**: Organize molecules into logical libraries based on project, target, or chemical class
2. **Use Consistent Naming**: Establish naming conventions for molecules and libraries
3. **Document Experiments**: Maintain clear records of experimental conditions and parameters
4. **Track Provenance**: Record the source and history of each molecule

## Quality Control

1. **Regular Validation**: Periodically validate molecular data for consistency
2. **Review Outliers**: Investigate molecules with extreme property values
3. **Update Properties**: Recalculate properties when new methods become available
4. **Verify Results**: Cross-check experimental results with predicted properties

# Troubleshooting

## Common Import Issues

| Issue | Possible Causes | Solutions |
|-------|----------------|----------|
| Invalid SMILES | Syntax errors, typos, unsupported features | Check SMILES syntax, verify with external tools, simplify complex structures |
| File Too Large | Exceeds 50 MB limit | Split into smaller files, remove unnecessary columns |
| Missing SMILES Column | Column not included or not mapped | Ensure CSV contains SMILES column, map correctly during import |
| Encoding Issues | Non-UTF-8 encoding | Save CSV with UTF-8 encoding |
| Performance Issues | Too many molecules or properties | Process in smaller batches, remove unnecessary columns |

## Validation Errors

| Error | Description | Resolution |
|-------|-------------|------------|
| Invalid Atom | Unrecognized atom symbol | Correct atom symbol, check for typos |
| Unbalanced Parentheses | Missing opening or closing parenthesis | Add missing parenthesis, verify structure |
| Unclosed Ring | Ring number used only once | Add matching ring number, verify structure |
| Invalid Valence | Atom has impossible number of bonds | Correct structure, check formal charges |
| Property Out of Range | Value exceeds expected range | Verify value, check calculation method |

## Performance Optimization

1. **Limit Columns**: Include only necessary properties in CSV files
2. **Pre-validate**: Clean data before import to reduce processing time
3. **Batch Processing**: Split very large datasets into multiple imports
4. **System Resources**: Ensure sufficient memory and CPU for large imports

# Appendix

## Property Calculation Methods

The platform uses RDKit // rdkit 2023.03+ for property calculations with the following specific methods:

| Property | RDKit Function | Notes |
|----------|---------------|-------|
| molecular_weight | Descriptors.MolWt | Includes all isotopes |
| logp | Descriptors.MolLogP | Crippen method |
| h_bond_donors | Descriptors.NumHDonors | Based on Lipinski definition |
| h_bond_acceptors | Descriptors.NumHAcceptors | Based on Lipinski definition |
| rotatable_bonds | Descriptors.NumRotatableBonds | Excludes terminal bonds |
| polar_surface_area | Descriptors.TPSA | Topological method |
| heavy_atom_count | Descriptors.HeavyAtomCount | Non-hydrogen atoms |
| ring_count | Descriptors.RingCount | SSSR definition |
| aromatic_rings | Lipinski.NumAromaticRings | RDKit aromaticity model |
| solubility | Custom calculation | Estimated from structure |

## SMILES Specification Reference

For detailed information about SMILES notation, refer to:

- Daylight SMILES specification: [http://www.daylight.com/smiles/](http://www.daylight.com/smiles/)
- OpenSMILES specification: [http://opensmiles.org/](http://opensmiles.org/)
- RDKit documentation: [https://www.rdkit.org/docs/](https://www.rdkit.org/docs/)

## Additional Resources

- Lipinski, C. A. (2004). Lead- and drug-like compounds: the rule-of-five revolution. Drug Discovery Today: Technologies, 1(4), 337-341.
- Veber, D. F., et al. (2002). Molecular properties that influence the oral bioavailability of drug candidates. Journal of Medicinal Chemistry, 45(12), 2615-2623.
- Wildman, S. A., & Crippen, G. M. (1999). Prediction of physicochemical parameters by atomic contributions. Journal of Chemical Information and Computer Sciences, 39(5), 868-873.