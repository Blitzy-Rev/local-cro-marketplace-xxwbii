"""
SQLAlchemy model imports for the Molecular Data Management and CRO Integration Platform.

This module centralizes imports of all SQLAlchemy models to ensure they are properly
registered with SQLAlchemy's metadata. It serves as a single import point for accessing
the Base class and all model definitions, which is particularly important for database
operations, migrations, and schema generation.

Importing this module ensures that all models are available to SQLAlchemy's metadata,
which is necessary for Alembic migrations to detect and generate proper schema changes.
"""

# Import the Base class
from .base_class import Base

# Import all models to register them with SQLAlchemy metadata
from ..models.user import User
from ..models.notification import Notification
from ..models.molecule import Molecule
from ..models.molecule_property import MoleculeProperty
from ..models.library import Library
from ..models.library_molecule import LibraryMolecule
from ..models.experiment_type import ExperimentType
from ..models.experiment_parameter import ExperimentParameter
from ..models.experiment import Experiment
from ..models.experiment_molecule import ExperimentMolecule
from ..models.submission import Submission
from ..models.submission_detail import SubmissionDetail
from ..models.result import Result
from ..models.result_file import ResultFile
from ..models.result_data import ResultData

# Export the Base class
__all__ = ["Base"]