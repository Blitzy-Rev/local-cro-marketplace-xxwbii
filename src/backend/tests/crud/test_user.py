import pytest
from datetime import datetime
import random
import string

# Internal imports
from ../../app/crud/crud_user import user
from ../../app/models/user import User
from ../../app/schemas/user import UserCreate, UserUpdate
from ../../app/core/security import verify_password
from ../../app/constants import UserRole, UserStatus


def test_create_user(db_session):
    """Test creating a new user with password"""
    email = f"test-{random.randint(1000, 9999)}@example.com"
    password = "SecureP@ssw0rd123"
    role = "pharma"
    
    # Create a UserCreate object with test data
    user_create = UserCreate(email=email, password=password, role=role)
    
    # Call user.create_with_password with db_session and UserCreate object
    created_user = user.create_with_password(db_session, user_create)
    
    # Assert that the returned user has the expected email and role
    assert created_user.email == email
    assert created_user.role == UserRole.PHARMA
    
    # Assert that the password was properly hashed (not stored in plaintext)
    assert created_user.password_hash != password
    assert verify_password(password, created_user.password_hash)
    
    # Assert that the user has default values set correctly
    assert created_user.is_active is True
    assert created_user.email_verified is False
    assert created_user.status == UserStatus.PENDING
    assert isinstance(created_user.created_at, datetime)
    assert created_user.last_login is None


def test_get_user_by_email(db_session, test_pharma_user):
    """Test retrieving a user by email"""
    # Get the email from the test_pharma_user fixture
    email = test_pharma_user.email
    
    # Call user.get_by_email with db_session and the email
    found_user = user.get_by_email(db_session, email)
    
    # Assert that the returned user is not None
    assert found_user is not None
    
    # Assert that the returned user has the expected email
    assert found_user.email == email
    
    # Test with a non-existent email and assert that None is returned
    non_existent_email = f"non-existent-{random.randint(1000, 9999)}@example.com"
    not_found_user = user.get_by_email(db_session, non_existent_email)
    assert not_found_user is None


def test_authenticate_user(db_session, test_password):
    """Test user authentication with email and password"""
    # Create a test user with a known password
    email = f"auth-test-{random.randint(1000, 9999)}@example.com"
    user_create = UserCreate(email=email, password=test_password, role="pharma")
    created_user = user.create_with_password(db_session, user_create)
    
    # Call user.authenticate with db_session, email, and correct password
    authenticated_user = user.authenticate(db_session, email=email, password=test_password)
    
    # Assert that the returned user is not None and has the expected email
    assert authenticated_user is not None
    assert authenticated_user.email == email
    
    # Call user.authenticate with db_session, email, and incorrect password
    incorrect_password = test_password + "wrong"
    not_authenticated = user.authenticate(db_session, email=email, password=incorrect_password)
    
    # Assert that None is returned for incorrect password
    assert not_authenticated is None
    
    # Call user.authenticate with db_session, non-existent email, and any password
    non_existent_email = f"non-existent-{random.randint(1000, 9999)}@example.com"
    not_found = user.authenticate(db_session, email=non_existent_email, password=test_password)
    
    # Assert that None is returned for non-existent user
    assert not_found is None


def test_update_user(db_session, test_pharma_user):
    """Test updating user information"""
    # Create a UserUpdate object with updated information
    new_email = f"updated-{random.randint(1000, 9999)}@example.com"
    user_update = UserUpdate(email=new_email, role="cro")
    
    # Call user.update_with_password with db_session, test_pharma_user, and UserUpdate object
    updated_user = user.update_with_password(db_session, test_pharma_user, user_update)
    
    # Assert that the returned user has the updated information
    assert updated_user.email == new_email
    assert updated_user.role == UserRole.CRO
    
    # Test updating with a dictionary instead of UserUpdate object
    update_dict = {"email": f"dict-update-{random.randint(1000, 9999)}@example.com", "status": "active"}
    dict_updated_user = user.update_with_password(db_session, updated_user, update_dict)
    
    # Assert that the returned user has the updated information
    assert dict_updated_user.email == update_dict["email"]
    assert dict_updated_user.status == UserStatus.ACTIVE


def test_update_password(db_session, test_pharma_user, test_password):
    """Test updating user password"""
    # Create a UserUpdate object with a new password
    new_password = f"NewP@ssw0rd{random.randint(1000, 9999)}"
    user_update = UserUpdate(password=new_password)
    
    # Store the current password hash for later comparison
    old_password_hash = test_pharma_user.password_hash
    
    # Call user.update_with_password with db_session, test_pharma_user, and UserUpdate object
    updated_user = user.update_with_password(db_session, test_pharma_user, user_update)
    
    # Assert that the password_hash has changed
    assert updated_user.password_hash != old_password_hash
    
    # Verify that the new password works with verify_password
    assert verify_password(new_password, updated_user.password_hash)
    
    # Authenticate with the new password and assert success
    authenticated = user.authenticate(db_session, email=updated_user.email, password=new_password)
    assert authenticated is not None
    assert authenticated.id == updated_user.id
    
    # Authenticate with the old password and assert failure
    not_authenticated = user.authenticate(db_session, email=updated_user.email, password=test_password)
    assert not_authenticated is None


def test_is_active(test_pharma_user):
    """Test checking if a user is active"""
    # Set test_pharma_user.is_active to True
    test_pharma_user.is_active = True
    
    # Call user.is_active with test_pharma_user
    assert user.is_active(test_pharma_user) is True
    
    # Set test_pharma_user.is_active to False
    test_pharma_user.is_active = False
    
    # Call user.is_active with test_pharma_user
    assert user.is_active(test_pharma_user) is False


def test_role_checks(test_pharma_user, test_cro_user, test_admin_user):
    """Test role checking functions"""
    # Call user.is_pharma with test_pharma_user and assert True
    assert user.is_pharma(test_pharma_user) is True
    # Call user.is_pharma with test_cro_user and assert False
    assert user.is_pharma(test_cro_user) is False
    # Call user.is_pharma with test_admin_user and assert False
    assert user.is_pharma(test_admin_user) is False
    
    # Call user.is_cro with test_cro_user and assert True
    assert user.is_cro(test_cro_user) is True
    # Call user.is_cro with test_pharma_user and assert False
    assert user.is_cro(test_pharma_user) is False
    # Call user.is_cro with test_admin_user and assert False
    assert user.is_cro(test_admin_user) is False
    
    # Call user.is_admin with test_admin_user and assert True
    assert user.is_admin(test_admin_user) is True
    # Call user.is_admin with test_pharma_user and assert False
    assert user.is_admin(test_pharma_user) is False
    # Call user.is_admin with test_cro_user and assert False
    assert user.is_admin(test_cro_user) is False


def test_update_last_login(db_session, test_pharma_user):
    """Test updating user's last login timestamp"""
    # Store the current last_login value (might be None)
    old_last_login = test_pharma_user.last_login
    
    # Call user.update_last_login with db_session and test_pharma_user
    updated_user = user.update_last_login(db_session, test_pharma_user)
    
    # Assert that last_login is now set and is a datetime object
    assert updated_user.last_login is not None
    assert isinstance(updated_user.last_login, datetime)
    
    # Assert that last_login is recent (within the last minute)
    assert (datetime.utcnow() - updated_user.last_login).total_seconds() < 60


def test_verify_email(db_session, test_pharma_user):
    """Test marking user's email as verified"""
    # Set test_pharma_user.email_verified to False
    test_pharma_user.email_verified = False
    db_session.add(test_pharma_user)
    db_session.commit()
    
    # Call user.verify_email with db_session and test_pharma_user
    verified_user = user.verify_email(db_session, test_pharma_user)
    
    # Assert that email_verified is now True
    assert verified_user.email_verified is True


def test_get_multi_by_role(db_session, test_pharma_user, test_cro_user, test_admin_user):
    """Test retrieving multiple users by role"""
    # Create additional pharma users for testing pagination
    additional_pharma_users = []
    for i in range(3):
        email = f"additional-pharma-{random.randint(1000, 9999)}@example.com"
        password = f"SecureP@ssw0rd{random.randint(1000, 9999)}"
        user_create = UserCreate(email=email, password=password, role="pharma")
        new_user = user.create_with_password(db_session, user_create)
        additional_pharma_users.append(new_user)
    
    # Call user.get_multi_by_role with db_session and 'pharma'
    pharma_users = user.get_multi_by_role(db_session, "pharma")
    
    # Assert that the returned list contains all pharma users
    assert len(pharma_users) >= 4  # test_pharma_user + 3 additional ones
    assert test_pharma_user.id in [u.id for u in pharma_users]
    for additional_user in additional_pharma_users:
        assert additional_user.id in [u.id for u in pharma_users]
    
    # Test pagination by setting skip and limit parameters
    limited_pharma_users = user.get_multi_by_role(db_session, "pharma", skip=1, limit=2)
    
    # Assert that the correct subset of users is returned
    assert len(limited_pharma_users) == 2
    
    # Test with other roles (cro, admin) and verify correct results
    cro_users = user.get_multi_by_role(db_session, "cro")
    assert len(cro_users) >= 1
    assert test_cro_user.id in [u.id for u in cro_users]
    
    admin_users = user.get_multi_by_role(db_session, "admin")
    assert len(admin_users) >= 1
    assert test_admin_user.id in [u.id for u in admin_users]


def test_count_by_role(db_session, test_pharma_user, test_cro_user, test_admin_user):
    """Test counting users by role"""
    # Create additional users with different roles
    initial_pharma_count = user.count_by_role(db_session, "pharma")
    initial_cro_count = user.count_by_role(db_session, "cro")
    initial_admin_count = user.count_by_role(db_session, "admin")
    
    # Create an additional user with each role
    for role in ["pharma", "cro", "admin"]:
        email = f"count-test-{role}-{random.randint(1000, 9999)}@example.com"
        password = f"SecureP@ssw0rd{random.randint(1000, 9999)}"
        user_create = UserCreate(email=email, password=password, role=role)
        user.create_with_password(db_session, user_create)
    
    # Call user.count_by_role with db_session and 'pharma'
    new_pharma_count = user.count_by_role(db_session, "pharma")
    new_cro_count = user.count_by_role(db_session, "cro")
    new_admin_count = user.count_by_role(db_session, "admin")
    
    # Assert that the count matches the expected number of pharma users
    assert new_pharma_count == initial_pharma_count + 1
    assert new_cro_count == initial_cro_count + 1
    assert new_admin_count == initial_admin_count + 1
    
    # Test with a role that has no users and verify count is 0
    non_existent_role_count = user.count_by_role(db_session, "non-existent-role")
    assert non_existent_role_count == 0


def test_inactive_user_authentication(db_session, test_pharma_user, test_password):
    """Test that inactive users cannot authenticate"""
    # Set test_pharma_user.is_active to False
    test_pharma_user.is_active = False
    db_session.add(test_pharma_user)
    db_session.commit()
    
    # Attempt to authenticate with correct credentials
    authenticated_user = user.authenticate(db_session, email=test_pharma_user.email, password=test_password)
    
    # Assert that None is returned despite correct credentials
    assert authenticated_user is None
    
    # Reset user to active state for other tests
    test_pharma_user.is_active = True
    db_session.add(test_pharma_user)
    db_session.commit()