from sqlalchemy import Column, String, Boolean, DateTime, Enum, Integer, Text, ARRAY
from sqlalchemy.orm import relationship, backref
from datetime import datetime

from ..db.base_class import Base
from ..constants import UserRole, UserStatus


class User(Base):
    """
    SQLAlchemy model representing a user account in the system with role-based access control.
    
    This model stores user authentication information, role, status, and maintains
    relationships to other entities in the system such as molecules, libraries,
    experiments, and notifications.
    
    Attributes:
        email (str): User's email address, must be unique
        password_hash (str): Hashed password for authentication
        role (UserRole): User's role (PHARMA, CRO, ADMIN)
        status (UserStatus): Account status (PENDING, ACTIVE, INACTIVE, LOCKED)
        is_active (bool): Whether the account is active
        email_verified (bool): Whether the email has been verified
        created_at (datetime): Account creation timestamp
        updated_at (datetime): Last update timestamp
        last_login (datetime): Last login timestamp
        password_history (list): Array of previous password hashes for password history validation
        molecules (relationship): Molecules created by this user
        libraries (relationship): Libraries created by this user
        experiments (relationship): Experiments created by this user
        cro_submissions (relationship): Submissions assigned to this CRO user
        notifications (relationship): Notifications for this user
    """
    # Basic user information
    email = Column(String(255), unique=True, index=True, nullable=False)
    password_hash = Column(String(255), nullable=False)
    role = Column(Enum(UserRole), nullable=False)
    status = Column(Enum(UserStatus), default=UserStatus.PENDING, nullable=False)
    is_active = Column(Boolean(), default=True)
    email_verified = Column(Boolean(), default=False)
    
    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)
    last_login = Column(DateTime, nullable=True)
    
    # Security - stores previous password hashes to prevent reuse (compliance requirement)
    password_history = Column(ARRAY(String), nullable=False)
    
    # Relationships
    molecules = relationship("Molecule", back_populates="created_by_user", 
                             cascade="all, delete-orphan")
    
    libraries = relationship("Library", back_populates="created_by_user", 
                            cascade="all, delete-orphan")
    
    experiments = relationship("Experiment", back_populates="created_by_user", 
                              cascade="all, delete-orphan")
    
    cro_submissions = relationship("Submission", back_populates="cro_user",
                                  foreign_keys="Submission.cro_id")
    
    notifications = relationship("Notification", back_populates="user", 
                                cascade="all, delete-orphan")