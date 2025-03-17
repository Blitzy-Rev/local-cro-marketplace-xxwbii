from sqlalchemy import Column, ForeignKey, String, Float, Index  # sqlalchemy version 2.0+
from sqlalchemy.orm import relationship  # sqlalchemy.orm version 2.0+
from typing import Optional
from uuid import UUID  # standard library

from ..db.base_class import Base
from .result import Result
from .molecule import Molecule


class ResultData(Base):
    """
    SQLAlchemy model representing a structured data point from an experimental result.
    
    This model stores individual data points or measurements from an experimental result,
    linking specific values to both a result and a molecule. It enables the storage and 
    retrieval of quantitative experimental data in a structured format.
    
    Attributes:
        id (int): Primary key (inherited from Base)
        result_id (int): Foreign key referencing the result this data belongs to
        molecule_id (int): Foreign key referencing the molecule this data is associated with
        data_name (str): Name or identifier of the data point (e.g., "IC50", "Solubility")
        data_value (float): Numerical value of the measurement
        data_unit (str): Unit of measurement (e.g., "nM", "mg/mL")
        
        result (Result): Relationship to the associated result
        molecule (Molecule): Relationship to the associated molecule
    """
    
    # Foreign keys
    result_id = Column(ForeignKey('results.id'), nullable=False, index=True)
    molecule_id = Column(ForeignKey('molecules.id'), nullable=False, index=True)
    
    # Data point details
    data_name = Column(String(100), nullable=False)
    data_value = Column(Float, nullable=True)
    data_unit = Column(String(50), nullable=True)
    
    # Relationships
    result = relationship("Result", back_populates="data_points")
    molecule = relationship("Molecule", back_populates="result_data")
    
    # Indexes for efficient querying
    __table_args__ = (
        Index('ix_result_data_result_molecule_name', 'result_id', 'molecule_id', 'data_name'),
    )