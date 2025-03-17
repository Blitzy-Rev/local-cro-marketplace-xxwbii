from .base import CRUDBase  # Import the base CRUD class for re-export
from .crud_user import CRUDUser, user  # Import the user CRUD class and singleton instance
from .crud_notification import CRUDNotification, notification  # Import the notification CRUD class and singleton instance
from .crud_molecule import CRUDMolecule, molecule  # Import the molecule CRUD class and singleton instance
from .crud_library import CRUDLibrary, library  # Import the library CRUD class and singleton instance
from .crud_experiment import CRUDExperiment, experiment  # Import the experiment CRUD class and singleton instance
from .crud_submission import CRUDSubmission, submission  # Import the submission CRUD class and singleton instance
from .crud_result import CRUDResult, result  # Import the result CRUD class and singleton instance

# Export all CRUD classes and instances
__all__ = [
    "CRUDBase",  # Base class for all CRUD operations
    "CRUDUser",  # CRUD class for user operations
    "user",  # Singleton instance of CRUDUser
    "CRUDNotification",  # CRUD class for notification operations
    "notification",  # Singleton instance of CRUDNotification
    "CRUDMolecule",  # CRUD class for molecule operations
    "molecule",  # Singleton instance of CRUDMolecule
    "CRUDLibrary",  # CRUD class for library operations
    "library",  # Singleton instance of CRUDLibrary
    "CRUDExperiment",  # CRUD class for experiment operations
    "experiment",  # Singleton instance of CRUDExperiment
    "CRUDSubmission",  # CRUD class for submission operations
    "submission",  # Singleton instance of CRUDSubmission
    "CRUDResult",  # CRUD class for result operations
    "result",  # Singleton instance of CRUDResult
]