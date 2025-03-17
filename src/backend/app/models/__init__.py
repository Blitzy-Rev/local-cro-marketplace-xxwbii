"""
Models package initialization module for the Molecular Data Management and CRO Integration Platform.

This module imports and re-exports all database model classes to make them accessible
through a single import point. This simplifies imports throughout the application
as multiple models can be imported with a single import statement.

The models in this package represent the core entities in the system's data model,
including users, molecules, libraries, experiments, submissions, and results.
These SQLAlchemy models provide the ORM mapping between Python objects and
database tables, defining both the schema and relationships.

Example:
    from app.models import User, Molecule, Library
    
    # Instead of multiple imports:
    # from app.models.user import User
    # from app.models.molecule import Molecule
    # from app.models.library import Library
"""

# Import all models for re-export
from .user import User
from .notification import Notification
from .molecule import Molecule
from .molecule_property import MoleculeProperty
from .library import Library
from .library_molecule import LibraryMolecule
from .experiment_type import ExperimentType
from .experiment_parameter import ExperimentParameter
from .experiment import Experiment
from .experiment_molecule import ExperimentMolecule
from .submission import Submission
from .submission_detail import SubmissionDetail
from .result import Result
from .result_file import ResultFile
from .result_data import ResultData

# Define the list of models that should be available when importing from this package
__all__ = [
    'User',
    'Notification',
    'Molecule',
    'MoleculeProperty',
    'Library',
    'LibraryMolecule',
    'ExperimentType',
    'ExperimentParameter',
    'Experiment',
    'ExperimentMolecule',
    'Submission',
    'SubmissionDetail',
    'Result',
    'ResultFile',
    'ResultData',
]