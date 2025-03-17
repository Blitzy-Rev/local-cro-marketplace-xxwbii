from sqlalchemy import Column, String, ForeignKey, Text, Index, Integer
from sqlalchemy.orm import relationship
from uuid import UUID  # standard library

from ..db.base_class import Base
from .submission import Submission


class SubmissionDetail(Base):
    """
    SQLAlchemy model representing additional details and parameters for a submission to a CRO.
    
    This model allows for storing flexible key-value pairs of metadata and configuration
    parameters associated with experiment submissions, such as specific assay conditions,
    concentration ranges, temperature settings, or any other submission-specific information.
    
    Attributes:
        id (int): Primary key (inherited from Base)
        submission_id (int): Foreign key referencing the parent submission
        detail_name (str): Name of the detail or parameter
        detail_value (str): Value of the detail or parameter
        
        submission (Submission): Relationship to the parent submission
    """
    
    # Foreign keys
    submission_id = Column(Integer, ForeignKey('submissions.id'), nullable=False)
    
    # Detail fields
    detail_name = Column(String(255), nullable=False)
    detail_value = Column(Text, nullable=True)
    
    # Relationships
    submission = relationship("Submission", back_populates="details")
    
    # Indexes for efficient lookups
    __table_args__ = (
        Index('idx_submission_detail_name', 'submission_id', 'detail_name'),
    )