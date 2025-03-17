from sqlalchemy import Column, ForeignKey, String, Float, Boolean, Index  # sqlalchemy version 2.0+
from sqlalchemy.orm import relationship  # sqlalchemy.orm version 2.0+
from typing import UUID

from ..db.base_class import Base


class MoleculeProperty(Base):
    """
    SQLAlchemy model representing a property of a molecule.
    
    This model provides a flexible schema to accommodate various types of
    molecular properties such as LogP, molecular weight, solubility, etc.
    Properties can be imported from CSV files or calculated by the system.
    
    The model supports efficient filtering and sorting of molecules based on
    their properties, which is essential for the molecule organization feature.
    
    Attributes:
        molecule_id (int): Foreign key reference to the molecules table
        property_name (str): Name of the property (e.g., "LogP", "MW", "Solubility")
        property_value (float): Numerical value of the property
        property_unit (str): Unit of measurement for the property (e.g., "g/mol", "mg/mL")
        is_calculated (bool): Indicates whether the property was calculated by the system or imported
        molecule (Molecule): Relationship back to the parent molecule
    """
    
    # Foreign key to the molecules table with cascade delete
    molecule_id = Column(Integer, ForeignKey('molecules.id', ondelete='CASCADE'), nullable=False, index=True)
    
    # Property identification and value
    property_name = Column(String(100), nullable=False)
    property_value = Column(Float, nullable=True)
    property_unit = Column(String(50), nullable=True)
    is_calculated = Column(Boolean, default=False, nullable=False)
    
    # Define relationship back to the molecule
    molecule = relationship("Molecule", back_populates="properties")
    
    # Create composite index for efficient property filtering
    __table_args__ = (
        Index('idx_molecule_property_name_value', molecule_id, property_name, property_value),
        Index('idx_property_name_value', property_name, property_value),
    )
    
    def __repr__(self):
        """
        String representation of the MoleculeProperty instance.
        
        Returns:
            str: String representation with property name, value and unit
        """
        return f"<MoleculeProperty(molecule_id={self.molecule_id}, " \
               f"property_name='{self.property_name}', " \
               f"property_value={self.property_value}, " \
               f"property_unit='{self.property_unit}')>"