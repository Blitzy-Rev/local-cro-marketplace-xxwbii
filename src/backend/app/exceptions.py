from fastapi import HTTPException, status  # version 0.95+

class BaseAppException(Exception):
    """
    Base exception class for all application-specific exceptions.
    
    This class provides a standard structure for application exceptions and
    includes functionality to convert to FastAPI HTTPException for API responses.
    """
    
    def __init__(self, message: str, status_code: int, details: dict = None):
        """
        Initialize the base exception with message, status code, and details.
        
        Args:
            message: Human-readable error message
            status_code: HTTP status code to return
            details: Additional error details for debugging
        """
        super().__init__(message)
        self.message = message
        self.status_code = status_code
        self.details = details or {}
        
    def to_http_exception(self) -> HTTPException:
        """
        Convert the application exception to a FastAPI HTTPException.
        
        Returns:
            HTTPException: FastAPI HTTP exception with appropriate status code and details
        """
        content = {
            "message": self.message,
            "details": self.details
        }
        return HTTPException(status_code=self.status_code, detail=content)


class AuthenticationException(BaseAppException):
    """
    Exception raised for authentication failures such as invalid credentials or token issues.
    """
    
    def __init__(self, message: str, details: dict = None):
        """
        Initialize authentication exception with message and details.
        
        Args:
            message: Human-readable error message
            details: Additional error details for debugging
        """
        super().__init__(message, status.HTTP_401_UNAUTHORIZED, details)


class AuthorizationException(BaseAppException):
    """
    Exception raised for authorization failures such as insufficient permissions.
    """
    
    def __init__(self, message: str, details: dict = None):
        """
        Initialize authorization exception with message and details.
        
        Args:
            message: Human-readable error message
            details: Additional error details for debugging
        """
        super().__init__(message, status.HTTP_403_FORBIDDEN, details)


class ValidationException(BaseAppException):
    """
    Exception raised for data validation failures such as invalid input formats or values.
    """
    
    def __init__(self, message: str, details: dict = None):
        """
        Initialize validation exception with message and details.
        
        Args:
            message: Human-readable error message
            details: Additional error details for debugging
        """
        super().__init__(message, status.HTTP_422_UNPROCESSABLE_ENTITY, details)


class ResourceNotFoundException(BaseAppException):
    """
    Exception raised when a requested resource is not found.
    """
    
    def __init__(self, message: str, details: dict = None):
        """
        Initialize resource not found exception with message and details.
        
        Args:
            message: Human-readable error message
            details: Additional error details for debugging
        """
        super().__init__(message, status.HTTP_404_NOT_FOUND, details)


class MolecularProcessingException(BaseAppException):
    """
    Exception raised for errors in molecular data processing such as invalid SMILES strings.
    """
    
    def __init__(self, message: str, details: dict = None):
        """
        Initialize molecular processing exception with message and details.
        
        Args:
            message: Human-readable error message
            details: Additional error details for debugging
        """
        super().__init__(message, status.HTTP_422_UNPROCESSABLE_ENTITY, details)


class CSVProcessingException(BaseAppException):
    """
    Exception raised for errors in CSV file processing such as invalid format or mapping.
    """
    
    def __init__(self, message: str, details: dict = None):
        """
        Initialize CSV processing exception with message and details.
        
        Args:
            message: Human-readable error message
            details: Additional error details for debugging
        """
        super().__init__(message, status.HTTP_422_UNPROCESSABLE_ENTITY, details)


class FileStorageException(BaseAppException):
    """
    Exception raised for errors in file storage operations such as upload or download failures.
    """
    
    def __init__(self, message: str, details: dict = None):
        """
        Initialize file storage exception with message and details.
        
        Args:
            message: Human-readable error message
            details: Additional error details for debugging
        """
        super().__init__(message, status.HTTP_500_INTERNAL_SERVER_ERROR, details)


class DatabaseException(BaseAppException):
    """
    Exception raised for database-related errors such as connection issues or constraint violations.
    """
    
    def __init__(self, message: str, details: dict = None):
        """
        Initialize database exception with message and details.
        
        Args:
            message: Human-readable error message
            details: Additional error details for debugging
        """
        super().__init__(message, status.HTTP_500_INTERNAL_SERVER_ERROR, details)


class ExperimentException(BaseAppException):
    """
    Exception raised for errors related to experiment management.
    """
    
    def __init__(self, message: str, details: dict = None):
        """
        Initialize experiment exception with message and details.
        
        Args:
            message: Human-readable error message
            details: Additional error details for debugging
        """
        super().__init__(message, status.HTTP_422_UNPROCESSABLE_ENTITY, details)


class SubmissionException(BaseAppException):
    """
    Exception raised for errors related to CRO submission processes.
    """
    
    def __init__(self, message: str, details: dict = None):
        """
        Initialize submission exception with message and details.
        
        Args:
            message: Human-readable error message
            details: Additional error details for debugging
        """
        super().__init__(message, status.HTTP_422_UNPROCESSABLE_ENTITY, details)


class ResultException(BaseAppException):
    """
    Exception raised for errors related to experimental results processing.
    """
    
    def __init__(self, message: str, details: dict = None):
        """
        Initialize result exception with message and details.
        
        Args:
            message: Human-readable error message
            details: Additional error details for debugging
        """
        super().__init__(message, status.HTTP_422_UNPROCESSABLE_ENTITY, details)


class ConfigurationException(BaseAppException):
    """
    Exception raised for system configuration errors.
    """
    
    def __init__(self, message: str, details: dict = None):
        """
        Initialize configuration exception with message and details.
        
        Args:
            message: Human-readable error message
            details: Additional error details for debugging
        """
        super().__init__(message, status.HTTP_500_INTERNAL_SERVER_ERROR, details)