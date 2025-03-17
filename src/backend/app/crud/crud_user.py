from datetime import datetime
from typing import Optional, List, Dict, Any, Union

from sqlalchemy.orm import Session  # sqlalchemy.orm version 2.0+
from sqlalchemy import select

from .base import CRUDBase
from ..models.user import User
from ..schemas.user import UserCreate, UserUpdate
from ..core.security import get_password_hash, verify_password
from ..constants import UserRole, UserStatus


class CRUDUser(CRUDBase[User, UserCreate, UserUpdate]):
    """
    CRUD operations for User model with specialized methods for user management
    """
    
    def __init__(self):
        """
        Initialize the CRUD user object with the User model
        """
        super().__init__(User)
    
    def get_by_email(self, db: Session, email: str) -> Optional[User]:
        """
        Get a user by email address
        
        Args:
            db: Database session
            email: Email address to search for
            
        Returns:
            User instance if found, None otherwise
        """
        stmt = select(User).where(User.email == email)
        result = db.execute(stmt).scalar_one_or_none()
        return result
    
    def create_with_password(self, db: Session, obj_in: UserCreate) -> User:
        """
        Create a new user with hashed password
        
        Args:
            db: Database session
            obj_in: User creation schema
            
        Returns:
            Created user instance
        """
        password_hash = get_password_hash(obj_in.password)
        
        # Map role string to enum
        role_map = {
            'pharma': UserRole.PHARMA,
            'cro': UserRole.CRO,
            'admin': UserRole.ADMIN
        }
        role = role_map.get(obj_in.role.lower(), UserRole.PHARMA)
        
        # Create user object with required fields
        db_obj_data = {
            "email": obj_in.email,
            "password_hash": password_hash,
            "role": role,
            "status": UserStatus.PENDING,
            "is_active": True,
            "email_verified": False,
            "created_at": datetime.utcnow(),
            "updated_at": datetime.utcnow(),
            "password_history": [password_hash]
        }
        
        db_obj = User(**db_obj_data)
        db.add(db_obj)
        db.commit()
        db.refresh(db_obj)
        return db_obj
    
    def update_with_password(
        self, db: Session, db_obj: User, obj_in: Union[UserUpdate, Dict[str, Any]]
    ) -> User:
        """
        Update user data, hashing password if provided
        
        Args:
            db: Database session
            db_obj: Existing user object
            obj_in: User update schema or dict with fields to update
            
        Returns:
            Updated user instance
        """
        update_data = obj_in if isinstance(obj_in, dict) else obj_in.dict(exclude_unset=True)
        
        if "password" in update_data:
            hashed_password = get_password_hash(update_data["password"])
            update_data["password_hash"] = hashed_password
            # Add to password history
            password_history = db_obj.password_history.copy() if db_obj.password_history else []
            update_data["password_history"] = [hashed_password] + password_history
            del update_data["password"]
        
        # Map role string to enum if provided
        if "role" in update_data:
            role_map = {
                'pharma': UserRole.PHARMA,
                'cro': UserRole.CRO,
                'admin': UserRole.ADMIN
            }
            update_data["role"] = role_map.get(update_data["role"].lower(), UserRole.PHARMA)
            
        # Map status string to enum if provided
        if "status" in update_data:
            status_map = {
                'pending': UserStatus.PENDING,
                'active': UserStatus.ACTIVE,
                'inactive': UserStatus.INACTIVE,
                'locked': UserStatus.LOCKED
            }
            update_data["status"] = status_map.get(update_data["status"].lower(), UserStatus.PENDING)
        
        update_data["updated_at"] = datetime.utcnow()
        
        return super().update(db, db_obj=db_obj, obj_in=update_data)
    
    def authenticate(self, db: Session, email: str, password: str) -> Optional[User]:
        """
        Authenticate a user by email and password
        
        Args:
            db: Database session
            email: User's email
            password: User's password
            
        Returns:
            Authenticated user if credentials are valid, None otherwise
        """
        user = self.get_by_email(db, email=email)
        if not user or not user.is_active:
            return None
        
        if not verify_password(password, user.password_hash):
            return None
        
        return user
    
    def is_active(self, user: User) -> bool:
        """
        Check if a user is active
        
        Args:
            user: User object
            
        Returns:
            True if user is active, False otherwise
        """
        return user.is_active
    
    def is_admin(self, user: User) -> bool:
        """
        Check if a user has admin role
        
        Args:
            user: User object
            
        Returns:
            True if user is admin, False otherwise
        """
        return user.role == UserRole.ADMIN
    
    def is_pharma(self, user: User) -> bool:
        """
        Check if a user has pharma role
        
        Args:
            user: User object
            
        Returns:
            True if user is pharma, False otherwise
        """
        return user.role == UserRole.PHARMA
    
    def is_cro(self, user: User) -> bool:
        """
        Check if a user has CRO role
        
        Args:
            user: User object
            
        Returns:
            True if user is CRO, False otherwise
        """
        return user.role == UserRole.CRO
    
    def update_last_login(self, db: Session, user: User) -> User:
        """
        Update user's last login timestamp
        
        Args:
            db: Database session
            user: User object
            
        Returns:
            Updated user instance
        """
        user.last_login = datetime.utcnow()
        db.add(user)
        db.commit()
        db.refresh(user)
        return user
    
    def verify_email(self, db: Session, user: User) -> User:
        """
        Mark user's email as verified
        
        Args:
            db: Database session
            user: User object
            
        Returns:
            Updated user instance
        """
        user.email_verified = True
        db.add(user)
        db.commit()
        db.refresh(user)
        return user
    
    def get_multi_by_role(
        self, db: Session, role: str, skip: int = 0, limit: int = 100
    ) -> List[User]:
        """
        Get multiple users filtered by role with pagination
        
        Args:
            db: Database session
            role: Role to filter by
            skip: Number of records to skip
            limit: Maximum number of records to return
            
        Returns:
            List of user instances with specified role
        """
        role_map = {
            'pharma': UserRole.PHARMA,
            'cro': UserRole.CRO,
            'admin': UserRole.ADMIN
        }
        role_enum = role_map.get(role.lower(), UserRole.PHARMA)
        
        stmt = select(User).where(User.role == role_enum).offset(skip).limit(limit)
        result = db.execute(stmt).scalars().all()
        return list(result)
    
    def count_by_role(self, db: Session, role: str) -> int:
        """
        Count users with a specific role
        
        Args:
            db: Database session
            role: Role to count
            
        Returns:
            Count of users with specified role
        """
        from sqlalchemy import func
        
        role_map = {
            'pharma': UserRole.PHARMA,
            'cro': UserRole.CRO,
            'admin': UserRole.ADMIN
        }
        role_enum = role_map.get(role.lower(), UserRole.PHARMA)
        
        stmt = select(func.count()).select_from(User).where(User.role == role_enum)
        result = db.execute(stmt).scalar_one()
        return result


# Create a singleton instance for application-wide use
user = CRUDUser()