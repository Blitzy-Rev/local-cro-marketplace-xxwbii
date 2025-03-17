"""
Unit tests for the CRUD operations related to CRO submissions in the Molecular Data Management
and CRO Integration Platform.

This file tests the functionality of the CRUDSubmission class, including submission creation,
retrieval, updating, status management, detail handling, and filtering operations.
"""

import pytest
from datetime import datetime, timedelta
import uuid

from app.crud.crud_submission import submission
from app.models.submission import Submission
from app.models.submission_detail import SubmissionDetail
from app.schemas.submission import SubmissionCreate, SubmissionUpdate
from app.constants import SubmissionStatus, ExperimentStatus


def test_create_submission(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test creating a new submission with details."""
    # Create a test experiment
    experiment = create_test_experiment(
        name="Test Experiment",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.DRAFT.name
    )
    
    # Create submission data
    submission_data = SubmissionCreate(
        experiment_id=experiment.id,
        cro_id=test_cro_user.id,
        details=[
            {"detail_name": "priority", "detail_value": "high"},
            {"detail_name": "notes", "detail_value": "Please process quickly"}
        ]
    )
    
    # Create submission
    new_submission = submission.create(db_session, submission_data)
    
    # Verify submission was created correctly
    assert new_submission.experiment_id == experiment.id
    assert new_submission.cro_id == test_cro_user.id
    assert new_submission.status == SubmissionStatus.PENDING.name
    assert new_submission.submitted_at is not None
    
    # Verify details were created
    details = submission.get_details(db_session, new_submission.id)
    assert len(details) == 2
    assert details["priority"] == "high"
    assert details["notes"] == "Please process quickly"


def test_get_submission_with_details(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test retrieving a submission with its details."""
    # Create a test experiment
    experiment = create_test_experiment(
        name="Test Experiment",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.DRAFT.name
    )
    
    # Create submission data
    submission_data = SubmissionCreate(
        experiment_id=experiment.id,
        cro_id=test_cro_user.id,
        details=[
            {"detail_name": "priority", "detail_value": "high"},
            {"detail_name": "notes", "detail_value": "Please process quickly"}
        ]
    )
    
    # Create submission
    new_submission = submission.create(db_session, submission_data)
    
    # Retrieve submission with details
    retrieved_submission = submission.get_with_details(db_session, new_submission.id)
    
    # Verify submission was retrieved correctly
    assert retrieved_submission.id == new_submission.id
    assert len(retrieved_submission.details) == 2
    
    # Verify details values
    detail_values = {detail.detail_name: detail.detail_value for detail in retrieved_submission.details}
    assert detail_values["priority"] == "high"
    assert detail_values["notes"] == "Please process quickly"


def test_get_submission_with_experiment(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test retrieving a submission with its associated experiment."""
    # Create a test experiment
    experiment = create_test_experiment(
        name="Test Experiment",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.DRAFT.name
    )
    
    # Create submission data
    submission_data = SubmissionCreate(
        experiment_id=experiment.id,
        cro_id=test_cro_user.id
    )
    
    # Create submission
    new_submission = submission.create(db_session, submission_data)
    
    # Retrieve submission with experiment
    retrieved_submission = submission.get_with_experiment(db_session, new_submission.id)
    
    # Verify submission was retrieved correctly
    assert retrieved_submission.id == new_submission.id
    assert retrieved_submission.experiment is not None
    assert retrieved_submission.experiment.id == experiment.id


def test_get_submission_with_cro(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test retrieving a submission with its associated CRO user."""
    # Create a test experiment
    experiment = create_test_experiment(
        name="Test Experiment",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.DRAFT.name
    )
    
    # Create submission data
    submission_data = SubmissionCreate(
        experiment_id=experiment.id,
        cro_id=test_cro_user.id
    )
    
    # Create submission
    new_submission = submission.create(db_session, submission_data)
    
    # Retrieve submission with CRO
    retrieved_submission = submission.get_with_cro(db_session, new_submission.id)
    
    # Verify submission was retrieved correctly
    assert retrieved_submission.id == new_submission.id
    assert retrieved_submission.cro is not None
    assert retrieved_submission.cro.id == test_cro_user.id


def test_update_submission(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test updating a submission's properties and details."""
    # Create a test experiment
    experiment = create_test_experiment(
        name="Test Experiment",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.DRAFT.name
    )
    
    # Create initial submission
    submission_data = SubmissionCreate(
        experiment_id=experiment.id,
        cro_id=test_cro_user.id,
        details=[
            {"detail_name": "priority", "detail_value": "normal"}
        ]
    )
    
    new_submission = submission.create(db_session, submission_data)
    
    # Update submission
    update_data = SubmissionUpdate(
        status=SubmissionStatus.QUOTE_PROVIDED.name,
        details=[
            {"detail_name": "priority", "detail_value": "high"},
            {"detail_name": "price", "detail_value": "1000"}
        ]
    )
    
    updated_submission = submission.update(db_session, new_submission, update_data)
    
    # Verify update was successful
    assert updated_submission.status == SubmissionStatus.QUOTE_PROVIDED.name
    assert updated_submission.updated_at > new_submission.submitted_at
    
    # Verify details were updated
    details = submission.get_details(db_session, updated_submission.id)
    assert len(details) == 2
    assert details["priority"] == "high"
    assert details["price"] == "1000"


def test_get_by_status(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test retrieving submissions by status."""
    # Create multiple test experiments
    experiment1 = create_test_experiment(
        name="Experiment 1",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.DRAFT.name
    )
    
    experiment2 = create_test_experiment(
        name="Experiment 2",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.DRAFT.name
    )
    
    # Create submissions with different statuses
    submission1 = submission.create(db_session, SubmissionCreate(
        experiment_id=experiment1.id,
        cro_id=test_cro_user.id,
        status=SubmissionStatus.PENDING.name
    ))
    
    submission2 = submission.create(db_session, SubmissionCreate(
        experiment_id=experiment2.id,
        cro_id=test_cro_user.id,
        status=SubmissionStatus.QUOTE_PROVIDED.name
    ))
    
    # Retrieve submissions by status
    pending_submissions = submission.get_by_status(db_session, SubmissionStatus.PENDING.name)
    quote_provided_submissions = submission.get_by_status(db_session, SubmissionStatus.QUOTE_PROVIDED.name)
    
    # Verify correct submissions were retrieved
    assert len(pending_submissions) == 1
    assert pending_submissions[0].id == submission1.id
    
    assert len(quote_provided_submissions) == 1
    assert quote_provided_submissions[0].id == submission2.id
    
    # Test count by status
    pending_count = submission.count_by_status(db_session, SubmissionStatus.PENDING.name)
    quote_provided_count = submission.count_by_status(db_session, SubmissionStatus.QUOTE_PROVIDED.name)
    
    assert pending_count == 1
    assert quote_provided_count == 1


def test_update_status(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test updating a submission's status."""
    # Create a test experiment
    experiment = create_test_experiment(
        name="Test Experiment",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.SUBMITTED.name
    )
    
    # Create submission
    new_submission = submission.create(db_session, SubmissionCreate(
        experiment_id=experiment.id,
        cro_id=test_cro_user.id,
        status=SubmissionStatus.PENDING.name
    ))
    
    # Update status
    updated_submission = submission.update_status(db_session, new_submission.id, SubmissionStatus.QUOTE_PROVIDED.name)
    
    # Verify status was updated
    assert updated_submission.status == SubmissionStatus.QUOTE_PROVIDED.name
    assert updated_submission.updated_at > new_submission.submitted_at
    
    # Verify experiment status was updated
    db_session.refresh(experiment)
    assert experiment.status == ExperimentStatus.QUOTE_PENDING.name


def test_get_by_experiment(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test retrieving submissions for a specific experiment."""
    # Create multiple test experiments
    experiment1 = create_test_experiment(
        name="Experiment 1",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.DRAFT.name
    )
    
    experiment2 = create_test_experiment(
        name="Experiment 2",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.DRAFT.name
    )
    
    # Create submissions for different experiments
    submission1 = submission.create(db_session, SubmissionCreate(
        experiment_id=experiment1.id,
        cro_id=test_cro_user.id
    ))
    
    submission2 = submission.create(db_session, SubmissionCreate(
        experiment_id=experiment2.id,
        cro_id=test_cro_user.id
    ))
    
    # Retrieve submissions by experiment
    exp1_submissions = submission.get_by_experiment(db_session, experiment1.id)
    exp2_submissions = submission.get_by_experiment(db_session, experiment2.id)
    
    # Verify correct submissions were retrieved
    assert len(exp1_submissions) == 1
    assert exp1_submissions[0].id == submission1.id
    
    assert len(exp2_submissions) == 1
    assert exp2_submissions[0].id == submission2.id
    
    # Test count by experiment
    exp1_count = submission.count_by_experiment(db_session, experiment1.id)
    exp2_count = submission.count_by_experiment(db_session, experiment2.id)
    
    assert exp1_count == 1
    assert exp2_count == 1


def test_get_by_cro(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test retrieving submissions assigned to a specific CRO."""
    # Create test experiments
    experiment1 = create_test_experiment(
        name="Experiment 1",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.DRAFT.name
    )
    
    experiment2 = create_test_experiment(
        name="Experiment 2",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.DRAFT.name
    )
    
    # Create a second CRO user (for testing purposes)
    from app.models.user import User
    from app.core.security import get_password_hash
    from app.constants import UserRole, UserStatus
    
    cro_user2 = User(
        email="cro2@example.com",
        password_hash=get_password_hash("Password123!"),
        role=UserRole.CRO,
        status=UserStatus.ACTIVE,
        is_active=True,
        email_verified=True,
        password_history=[]
    )
    db_session.add(cro_user2)
    db_session.flush()
    
    # Create submissions for different CROs
    submission1 = submission.create(db_session, SubmissionCreate(
        experiment_id=experiment1.id,
        cro_id=test_cro_user.id
    ))
    
    submission2 = submission.create(db_session, SubmissionCreate(
        experiment_id=experiment2.id,
        cro_id=cro_user2.id
    ))
    
    # Retrieve submissions by CRO
    cro1_submissions = submission.get_by_cro(db_session, test_cro_user.id)
    cro2_submissions = submission.get_by_cro(db_session, cro_user2.id)
    
    # Verify correct submissions were retrieved
    assert len(cro1_submissions) == 1
    assert cro1_submissions[0].id == submission1.id
    
    assert len(cro2_submissions) == 1
    assert cro2_submissions[0].id == submission2.id
    
    # Test count by CRO
    cro1_count = submission.count_by_cro(db_session, test_cro_user.id)
    cro2_count = submission.count_by_cro(db_session, cro_user2.id)
    
    assert cro1_count == 1
    assert cro2_count == 1


def test_submission_detail_operations(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test adding, retrieving, and removing submission details."""
    # Create a test experiment
    experiment = create_test_experiment(
        name="Test Experiment",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.DRAFT.name
    )
    
    # Create submission without details
    new_submission = submission.create(db_session, SubmissionCreate(
        experiment_id=experiment.id,
        cro_id=test_cro_user.id
    ))
    
    # Add detail
    detail = submission.add_detail(db_session, new_submission.id, "priority", "high")
    
    # Verify detail was added
    assert detail.submission_id == new_submission.id
    assert detail.detail_name == "priority"
    assert detail.detail_value == "high"
    
    # Get detail
    priority_value = submission.get_detail(db_session, new_submission.id, "priority")
    assert priority_value == "high"
    
    # Get all details
    all_details = submission.get_details(db_session, new_submission.id)
    assert "priority" in all_details
    assert all_details["priority"] == "high"
    
    # Remove detail
    result = submission.remove_detail(db_session, new_submission.id, "priority")
    assert result is True
    
    # Verify detail was removed
    priority_value = submission.get_detail(db_session, new_submission.id, "priority")
    assert priority_value is None


def test_provide_quote(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test providing a quote for a submission."""
    # Create a test experiment
    experiment = create_test_experiment(
        name="Test Experiment",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.SUBMITTED.name
    )
    
    # Create submission
    new_submission = submission.create(db_session, SubmissionCreate(
        experiment_id=experiment.id,
        cro_id=test_cro_user.id,
        status=SubmissionStatus.PENDING.name
    ))
    
    # Provide quote
    updated_submission = submission.provide_quote(
        db_session, 
        new_submission.id, 
        price=1000.00, 
        turnaround_time=7, 
        notes="Standard pricing for this type of experiment."
    )
    
    # Verify submission was updated
    assert updated_submission.status == SubmissionStatus.QUOTE_PROVIDED.name
    assert updated_submission.price == 1000.00
    assert updated_submission.turnaround_days == 7
    
    # Verify details were added
    details = submission.get_details(db_session, updated_submission.id)
    assert details["price"] == "1000.0"
    assert details["currency"] == "USD"
    assert details["turnaround_days"] == "7"
    assert details["quote_notes"] == "Standard pricing for this type of experiment."
    
    # Verify experiment status was updated
    db_session.refresh(experiment)
    assert experiment.status == ExperimentStatus.QUOTE_PENDING.name


def test_respond_to_quote_approved(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test approving a quote for a submission."""
    # Create a test experiment
    experiment = create_test_experiment(
        name="Test Experiment",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.QUOTE_PENDING.name
    )
    
    # Create submission and provide quote
    new_submission = submission.create(db_session, SubmissionCreate(
        experiment_id=experiment.id,
        cro_id=test_cro_user.id,
        status=SubmissionStatus.PENDING.name
    ))
    
    quoted_submission = submission.provide_quote(
        db_session, 
        new_submission.id, 
        price=1000.00, 
        turnaround_time=7
    )
    
    # Approve quote
    updated_submission = submission.respond_to_quote(
        db_session,
        quoted_submission.id,
        approved=True,
        notes="Quote approved, please proceed."
    )
    
    # Verify submission was updated
    assert updated_submission.status == SubmissionStatus.APPROVED.name
    
    # Verify response notes were added
    response_notes = submission.get_detail(db_session, updated_submission.id, "response_notes")
    assert response_notes == "Quote approved, please proceed."
    
    # Verify experiment status was updated
    db_session.refresh(experiment)
    assert experiment.status == ExperimentStatus.IN_PROGRESS.name


def test_respond_to_quote_rejected(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test rejecting a quote for a submission."""
    # Create a test experiment
    experiment = create_test_experiment(
        name="Test Experiment",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.QUOTE_PENDING.name
    )
    
    # Create submission and provide quote
    new_submission = submission.create(db_session, SubmissionCreate(
        experiment_id=experiment.id,
        cro_id=test_cro_user.id,
        status=SubmissionStatus.PENDING.name
    ))
    
    quoted_submission = submission.provide_quote(
        db_session, 
        new_submission.id, 
        price=1000.00, 
        turnaround_time=7
    )
    
    # Reject quote
    updated_submission = submission.respond_to_quote(
        db_session,
        quoted_submission.id,
        approved=False,
        notes="Price is too high."
    )
    
    # Verify submission was updated
    assert updated_submission.status == SubmissionStatus.QUOTE_REJECTED.name
    
    # Verify response notes were added
    response_notes = submission.get_detail(db_session, updated_submission.id, "response_notes")
    assert response_notes == "Price is too high."
    
    # Verify experiment status was updated
    db_session.refresh(experiment)
    assert experiment.status == ExperimentStatus.QUOTE_REJECTED.name


def test_filter_submissions(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test filtering submissions based on various criteria."""
    # Create multiple test experiments
    experiment1 = create_test_experiment(
        name="Experiment 1",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.DRAFT.name
    )
    
    experiment2 = create_test_experiment(
        name="Experiment 2",
        created_by=test_pharma_user.id,
        experiment_type_id=1,
        status=ExperimentStatus.DRAFT.name
    )
    
    # Create submissions with various properties
    submission1 = submission.create(db_session, SubmissionCreate(
        experiment_id=experiment1.id,
        cro_id=test_cro_user.id,
        status=SubmissionStatus.PENDING.name
    ))
    
    # Manually update submitted_at to create date variation for filtering tests
    now = datetime.utcnow()
    yesterday = now - timedelta(days=1)
    submission1.submitted_at = yesterday
    db_session.add(submission1)
    db_session.commit()
    
    submission2 = submission.create(db_session, SubmissionCreate(
        experiment_id=experiment2.id,
        cro_id=test_cro_user.id,
        status=SubmissionStatus.QUOTE_PROVIDED.name
    ))
    
    # Create a second CRO user for testing
    from app.models.user import User
    from app.core.security import get_password_hash
    from app.constants import UserRole, UserStatus
    
    cro_user2 = User(
        email="cro2_filter@example.com",  # Different email from test_get_by_cro
        password_hash=get_password_hash("Password123!"),
        role=UserRole.CRO,
        status=UserStatus.ACTIVE,
        is_active=True,
        email_verified=True,
        password_history=[]
    )
    db_session.add(cro_user2)
    db_session.flush()
    
    submission3 = submission.create(db_session, SubmissionCreate(
        experiment_id=experiment1.id,
        cro_id=cro_user2.id,
        status=SubmissionStatus.PENDING.name
    ))
    
    # Test filtering by status
    status_filter = {"status": SubmissionStatus.PENDING.name}
    pending_submissions, count = submission.filter_submissions(db_session, status_filter)
    assert count == 2
    assert len(pending_submissions) == 2
    
    # Test filtering by experiment_id
    experiment_filter = {"experiment_id": experiment1.id}
    exp1_submissions, count = submission.filter_submissions(db_session, experiment_filter)
    assert count == 2
    assert len(exp1_submissions) == 2
    
    # Test filtering by cro_id
    cro_filter = {"cro_id": test_cro_user.id}
    cro1_submissions, count = submission.filter_submissions(db_session, cro_filter)
    assert count == 2
    assert len(cro1_submissions) == 2
    
    # Test filtering by date range
    date_filter = {"submitted_after": yesterday - timedelta(hours=1), "submitted_before": now - timedelta(hours=1)}
    date_submissions, count = submission.filter_submissions(db_session, date_filter)
    assert count >= 1  # At least one submission should be in this range
    
    # Test multiple filters
    combined_filter = {
        "experiment_id": experiment1.id,
        "cro_id": test_cro_user.id,
        "status": SubmissionStatus.PENDING.name
    }
    filtered_submissions, count = submission.filter_submissions(db_session, combined_filter)
    assert count == 1
    assert len(filtered_submissions) == 1
    assert filtered_submissions[0].id == submission1.id