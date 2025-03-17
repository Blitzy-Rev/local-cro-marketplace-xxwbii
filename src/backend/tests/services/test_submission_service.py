"""
Unit tests for the submission service in the Molecular Data Management and CRO Integration Platform.
This file contains comprehensive test cases for all submission service functions, including
submission creation, status management, quote handling, and CRO interactions.
"""

import pytest  # pytest version: latest
from pytest import fixture, mark, raises  # pytest version: latest
import unittest.mock  # unittest.mock version: standard library
from unittest.mock import MagicMock, patch  # unittest.mock version: standard library
import typing  # typing version: standard library
from typing import Dict, List, Any, Optional, UUID  # typing version: standard library
import uuid  # uuid version: standard library
from uuid import uuid4  # uuid version: standard library
import datetime  # datetime version: standard library
from datetime import datetime  # datetime version: standard library

from ..conftest import db_session  # Path: src/backend/tests/conftest.py
from ..conftest import test_pharma_user, test_cro_user  # Path: src/backend/tests/conftest.py
from ..conftest import create_test_experiment  # Path: src/backend/tests/conftest.py
from ...app.services.submission_service import SubmissionService  # Path: src/backend/app/services/submission_service.py
from ...app.services.submission_service import submission_service  # Path: src/backend/app/services/submission_service.py
from ...app.services.experiment_service import experiment_service  # Path: src/backend/app/services/experiment_service.py
from ...app.schemas.submission import SubmissionCreate  # Path: src/backend/app/schemas/submission.py
from ...app.schemas.submission import SubmissionUpdate  # Path: src/backend/app/schemas/submission.py
from ...app.schemas.submission import QuoteProvide  # Path: src/backend/app/schemas/submission.py
from ...app.schemas.submission import QuoteResponse  # Path: src/backend/app/schemas/submission.py
from ...app.schemas.submission import SubmissionFilter  # Path: src/backend/app/schemas/submission.py
from ...app.constants import SubmissionStatus  # Path: src/backend/app/constants.py
from ...app.exceptions import ValidationException, ResourceNotFoundException, SubmissionException  # Path: src/backend/app/exceptions.py


def test_create_submission_success(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test successful creation of a submission"""
    # Create a test experiment using create_test_experiment fixture
    experiment = create_test_experiment(name="Test Experiment", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    
    # Mock experiment_service.prepare_for_submission to return experiment data
    with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment.id}):
        # Create SubmissionCreate object with experiment_id and cro_id
        submission_create = SubmissionCreate(experiment_id=experiment.id, cro_id=test_cro_user.id)
        
        # Call submission_service.create_submission with the data and pharma user ID
        submission = submission_service.create_submission(submission_create, test_pharma_user.id)
        
        # Assert that the returned submission has the correct experiment_id and cro_id
        assert submission["experiment_id"] == experiment.id
        assert submission["cro_id"] == test_cro_user.id
        
        # Assert that the status is set to PENDING
        assert submission["status"] == SubmissionStatus.PENDING.name
        
        # Assert that submitted_at is set to a datetime
        assert isinstance(submission["submitted_at"], str)


def test_create_submission_experiment_not_ready(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test submission creation fails when experiment is not ready"""
    # Create a test experiment using create_test_experiment fixture
    experiment = create_test_experiment(name="Test Experiment", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    
    # Mock experiment_service.prepare_for_submission to raise ValidationException
    with patch.object(experiment_service, "get_experiment_by_id", side_effect=ResourceNotFoundException("Experiment not found")):
        # Create SubmissionCreate object with experiment_id and cro_id
        submission_create = SubmissionCreate(experiment_id=experiment.id, cro_id=test_cro_user.id)
        
        # Use pytest.raises to assert that ValidationException is raised
        with pytest.raises(ResourceNotFoundException) as exc_info:
            # Call submission_service.create_submission with the data and pharma user ID
            submission_service.create_submission(submission_create, test_pharma_user.id)
        
        # Assert that the exception has the expected message
        assert str(exc_info.value) == "Experiment not found"


def test_update_submission_success(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test successful update of a submission"""
    # Create a test experiment using create_test_experiment fixture
    experiment = create_test_experiment(name="Test Experiment", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    
    # Mock experiment_service.prepare_for_submission to return experiment data
    with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment.id}):
        # Create a submission using submission_service.create_submission
        submission_create = SubmissionCreate(experiment_id=experiment.id, cro_id=test_cro_user.id)
        submission = submission_service.create_submission(submission_create, test_pharma_user.id)
        
        # Create SubmissionUpdate object with updated notes
        submission_update = SubmissionUpdate(notes="Updated notes")
        
        # Call submission_service.update_submission with the submission ID, update data, and pharma user ID
        updated_submission = submission_service.update_submission(submission["id"], submission_update, test_pharma_user.id)
        
        # Assert that the returned submission has the updated notes
        assert updated_submission["notes"] == "Updated notes"
        
        # Assert that updated_at is set to a datetime
        assert isinstance(updated_submission["updated_at"], str)


def test_update_submission_not_found(db_session, test_pharma_user):
    """Test submission update fails when submission is not found"""
    # Generate a random UUID for a non-existent submission
    submission_id = uuid4()
    
    # Create SubmissionUpdate object with updated notes
    submission_update = SubmissionUpdate(notes="Updated notes")
    
    # Use pytest.raises to assert that ResourceNotFoundException is raised
    with pytest.raises(ResourceNotFoundException) as exc_info:
        # Call submission_service.update_submission with the random UUID, update data, and pharma user ID
        submission_service.update_submission(submission_id, submission_update, test_pharma_user.id)
    
    # Assert that the exception has the expected message
    assert str(exc_info.value) == f"Submission not found: {submission_id}"


def test_cancel_submission_success(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test successful cancellation of a submission"""
    # Create a test experiment using create_test_experiment fixture
    experiment = create_test_experiment(name="Test Experiment", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    
    # Mock experiment_service.prepare_for_submission to return experiment data
    with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment.id}):
        # Create a submission using submission_service.create_submission
        submission_create = SubmissionCreate(experiment_id=experiment.id, cro_id=test_cro_user.id)
        submission = submission_service.create_submission(submission_create, test_pharma_user.id)
        
        # Call submission_service.cancel_submission with the submission ID and pharma user ID
        submission_service.cancel_submission(submission["id"], test_pharma_user.id)
        
        # Get the updated submission using submission_service.get_submission_by_id
        updated_submission = submission_service.get_submission_by_id(submission["id"])
        
        # Assert that the submission status is CANCELLED
        assert updated_submission["status"] == SubmissionStatus.CANCELLED.name


def test_cancel_submission_not_found(db_session, test_pharma_user):
    """Test submission cancellation fails when submission is not found"""
    # Generate a random UUID for a non-existent submission
    submission_id = uuid4()
    
    # Use pytest.raises to assert that ResourceNotFoundException is raised
    with pytest.raises(ResourceNotFoundException) as exc_info:
        # Call submission_service.cancel_submission with the random UUID and pharma user ID
        submission_service.cancel_submission(submission_id, test_pharma_user.id)
    
    # Assert that the exception has the expected message
    assert str(exc_info.value) == f"Submission not found: {submission_id}"


def test_update_submission_status_success(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test successful update of submission status"""
    # Create a test experiment using create_test_experiment fixture
    experiment = create_test_experiment(name="Test Experiment", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    
    # Mock experiment_service.prepare_for_submission to return experiment data
    with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment.id}):
        # Create a submission using submission_service.create_submission
        submission_create = SubmissionCreate(experiment_id=experiment.id, cro_id=test_cro_user.id)
        submission = submission_service.create_submission(submission_create, test_pharma_user.id)
        
        # Call submission_service.update_submission_status with the submission ID, new status, and pharma user ID
        new_status = SubmissionStatus.QUOTE_PROVIDED.name
        submission_service.update_submission_status(submission["id"], new_status, test_pharma_user.id)
        
        # Get the updated submission using submission_service.get_submission_by_id
        updated_submission = submission_service.get_submission_by_id(submission["id"])
        
        # Assert that the submission status is updated to the new status
        assert updated_submission["status"] == new_status


def test_provide_quote_success(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test successful provision of a quote by CRO"""
    # Create a test experiment using create_test_experiment fixture
    experiment = create_test_experiment(name="Test Experiment", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    
    # Mock experiment_service.prepare_for_submission to return experiment data
    with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment.id}):
        # Create a submission using submission_service.create_submission
        submission_create = SubmissionCreate(experiment_id=experiment.id, cro_id=test_cro_user.id)
        submission = submission_service.create_submission(submission_create, test_pharma_user.id)
        
        # Create QuoteProvide object with price and turnaround_days
        quote_data = QuoteProvide(price=1000.00, turnaround_days=7)
        
        # Call submission_service.provide_quote with the submission ID, quote data, and CRO user ID
        submission_service.provide_quote(submission["id"], quote_data, test_cro_user.id)
        
        # Get the updated submission using submission_service.get_submission_by_id
        updated_submission = submission_service.get_submission_by_id(submission["id"])
        
        # Assert that the submission status is QUOTE_PROVIDED
        assert updated_submission["status"] == SubmissionStatus.QUOTE_PROVIDED.name
        
        # Assert that the price and turnaround_days are set correctly
        assert updated_submission["price"] == 1000.00
        assert updated_submission["turnaround_days"] == 7


def test_provide_quote_not_cro(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test quote provision fails when user is not the assigned CRO"""
    # Create a test experiment using create_test_experiment fixture
    experiment = create_test_experiment(name="Test Experiment", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    
    # Mock experiment_service.prepare_for_submission to return experiment data
    with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment.id}):
        # Create a submission using submission_service.create_submission
        submission_create = SubmissionCreate(experiment_id=experiment.id, cro_id=test_cro_user.id)
        submission = submission_service.create_submission(submission_create, test_pharma_user.id)
        
        # Create QuoteProvide object with price and turnaround_days
        quote_data = QuoteProvide(price=1000.00, turnaround_days=7)
        
        # Use pytest.raises to assert that AuthorizationException is raised
        with pytest.raises(Exception):
            # Call submission_service.provide_quote with the submission ID, quote data, and pharma user ID (not CRO)
            submission_service.provide_quote(submission["id"], quote_data, test_pharma_user.id)


def test_respond_to_quote_approve_success(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test successful approval of a quote by pharma user"""
    # Create a test experiment using create_test_experiment fixture
    experiment = create_test_experiment(name="Test Experiment", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    
    # Mock experiment_service.prepare_for_submission to return experiment data
    with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment.id}):
        # Create a submission using submission_service.create_submission
        submission_create = SubmissionCreate(experiment_id=experiment.id, cro_id=test_cro_user.id)
        submission = submission_service.create_submission(submission_create, test_pharma_user.id)
        
        # Update submission status to QUOTE_PROVIDED
        with patch.object(submission_service, "get_submission_by_id", return_value={"id": submission["id"], "status": SubmissionStatus.QUOTE_PROVIDED.name}):
            # Create QuoteResponse object with approved=True
            response_data = QuoteResponse(approved=True)
            
            # Call submission_service.respond_to_quote with the submission ID, response data, and pharma user ID
            submission_service.respond_to_quote(submission["id"], response_data, test_pharma_user.id)
            
            # Get the updated submission using submission_service.get_submission_by_id
            updated_submission = submission_service.get_submission_by_id(submission["id"])
            
            # Assert that the submission status is APPROVED
            assert updated_submission["status"] == SubmissionStatus.APPROVED.name


def test_respond_to_quote_reject_success(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test successful rejection of a quote by pharma user"""
    # Create a test experiment using create_test_experiment fixture
    experiment = create_test_experiment(name="Test Experiment", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    
    # Mock experiment_service.prepare_for_submission to return experiment data
    with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment.id}):
        # Create a submission using submission_service.create_submission
        submission_create = SubmissionCreate(experiment_id=experiment.id, cro_id=test_cro_user.id)
        submission = submission_service.create_submission(submission_create, test_pharma_user.id)
        
        # Update submission status to QUOTE_PROVIDED
        with patch.object(submission_service, "get_submission_by_id", return_value={"id": submission["id"], "status": SubmissionStatus.QUOTE_PROVIDED.name}):
            # Create QuoteResponse object with approved=False
            response_data = QuoteResponse(approved=False)
            
            # Call submission_service.respond_to_quote with the submission ID, response data, and pharma user ID
            submission_service.respond_to_quote(submission["id"], response_data, test_pharma_user.id)
            
            # Get the updated submission using submission_service.get_submission_by_id
            updated_submission = submission_service.get_submission_by_id(submission["id"])
            
            # Assert that the submission status is QUOTE_REJECTED
            assert updated_submission["status"] == SubmissionStatus.QUOTE_REJECTED.name


def test_get_submissions_filter_by_status(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test filtering submissions by status"""
    # Create multiple test experiments using create_test_experiment fixture
    experiment1 = create_test_experiment(name="Test Experiment 1", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    experiment2 = create_test_experiment(name="Test Experiment 2", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    
    # Mock experiment_service.prepare_for_submission to return experiment data
    with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment1.id}):
        # Create multiple submissions with different statuses
        submission_create1 = SubmissionCreate(experiment_id=experiment1.id, cro_id=test_cro_user.id, status=SubmissionStatus.PENDING.name)
        submission1 = submission_service.create_submission(submission_create1, test_pharma_user.id)
        
        with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment2.id}):
            submission_create2 = SubmissionCreate(experiment_id=experiment2.id, cro_id=test_cro_user.id, status=SubmissionStatus.QUOTE_PROVIDED.name)
            submission2 = submission_service.create_submission(submission_create2, test_pharma_user.id)
        
        # Create SubmissionFilter object with status=PENDING
        submission_filter = SubmissionFilter(status=SubmissionStatus.PENDING.name)
        
        # Call submission_service.get_submissions with the filter
        submissions = submission_service.get_submissions(submission_filter)
        
        # Assert that only submissions with PENDING status are returned
        assert len(submissions["items"]) == 1
        assert submissions["items"][0]["status"] == SubmissionStatus.PENDING.name
        
        # Assert that the total count matches the expected number
        assert submissions["total"] == 1


def test_get_submissions_by_experiment(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test retrieving submissions for a specific experiment"""
    # Create multiple test experiments using create_test_experiment fixture
    experiment1 = create_test_experiment(name="Test Experiment 1", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    experiment2 = create_test_experiment(name="Test Experiment 2", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    
    # Mock experiment_service.prepare_for_submission to return experiment data
    with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment1.id}):
        # Create submissions for different experiments
        submission_create1 = SubmissionCreate(experiment_id=experiment1.id, cro_id=test_cro_user.id)
        submission1 = submission_service.create_submission(submission_create1, test_pharma_user.id)
        
        with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment2.id}):
            submission_create2 = SubmissionCreate(experiment_id=experiment2.id, cro_id=test_cro_user.id)
            submission2 = submission_service.create_submission(submission_create2, test_pharma_user.id)
        
        # Call submission_service.get_submissions_by_experiment with the first experiment ID and pharma user ID
        submissions = submission_service.get_submissions_by_experiment(experiment1.id, test_pharma_user.id)
        
        # Assert that only submissions for the specified experiment are returned
        assert len(submissions["items"]) == 1
        assert submissions["items"][0]["experiment_id"] == experiment1.id
        
        # Assert that the total count matches the expected number
        assert submissions["total"] == 1


def test_get_submissions_by_cro(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test retrieving submissions assigned to a specific CRO"""
    # Create multiple test experiments using create_test_experiment fixture
    experiment1 = create_test_experiment(name="Test Experiment 1", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    experiment2 = create_test_experiment(name="Test Experiment 2", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    
    # Mock experiment_service.prepare_for_submission to return experiment data
    with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment1.id}):
        # Create submissions assigned to different CROs
        submission_create1 = SubmissionCreate(experiment_id=experiment1.id, cro_id=test_cro_user.id)
        submission1 = submission_service.create_submission(submission_create1, test_pharma_user.id)
        
        with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment2.id}):
            submission_create2 = SubmissionCreate(experiment_id=experiment2.id, cro_id=uuid4())
            submission2 = submission_service.create_submission(submission_create2, test_pharma_user.id)
        
        # Call submission_service.get_submissions_by_cro with the CRO user ID
        submissions = submission_service.get_submissions_by_cro(test_cro_user.id)
        
        # Assert that only submissions assigned to the specified CRO are returned
        assert len(submissions["items"]) == 1
        assert submissions["items"][0]["cro_id"] == test_cro_user.id
        
        # Assert that the total count matches the expected number
        assert submissions["total"] == 1


def test_get_submissions_by_status(db_session, test_pharma_user, test_cro_user, create_test_experiment):
    """Test retrieving submissions with a specific status"""
    # Create multiple test experiments using create_test_experiment fixture
    experiment1 = create_test_experiment(name="Test Experiment 1", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    experiment2 = create_test_experiment(name="Test Experiment 2", created_by=test_pharma_user.id, experiment_type_id=uuid4())
    
    # Mock experiment_service.prepare_for_submission to return experiment data
    with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment1.id}):
        # Create submissions with different statuses
        submission_create1 = SubmissionCreate(experiment_id=experiment1.id, cro_id=test_cro_user.id, status=SubmissionStatus.PENDING.name)
        submission1 = submission_service.create_submission(submission_create1, test_pharma_user.id)
        
        with patch.object(experiment_service, "get_experiment_by_id", return_value={"id": experiment2.id}):
            submission_create2 = SubmissionCreate(experiment_id=experiment2.id, cro_id=test_cro_user.id, status=SubmissionStatus.QUOTE_PROVIDED.name)
            submission2 = submission_service.create_submission(submission_create2, test_pharma_user.id)
        
        # Call submission_service.get_submissions_by_status with PENDING status
        submissions = submission_service.get_submissions_by_status(SubmissionStatus.PENDING.name)
        
        # Assert that only submissions with PENDING status are returned
        assert len(submissions["items"]) == 1
        assert submissions["items"][0]["status"] == SubmissionStatus.PENDING.name
        
        # Assert that the total count matches the expected number
        assert submissions["total"] == 1


def test_validate_status_transition_valid():
    """Test validation of valid status transitions"""
    # Define valid status transitions to test
    valid_transitions = [
        (SubmissionStatus.PENDING.name, SubmissionStatus.QUOTE_PROVIDED.name),
        (SubmissionStatus.QUOTE_PROVIDED.name, SubmissionStatus.APPROVED.name),
        (SubmissionStatus.APPROVED.name, SubmissionStatus.IN_PROGRESS.name),
        (SubmissionStatus.IN_PROGRESS.name, SubmissionStatus.COMPLETED.name)
    ]
    
    # For each transition, call submission_service.validate_status_transition
    for current_status, new_status in valid_transitions:
        # Assert that the result is True for each valid transition
        assert submission_service.validate_status_transition(current_status, new_status) is None


def test_validate_status_transition_invalid():
    """Test validation of invalid status transitions"""
    # Define invalid status transitions to test
    invalid_transitions = [
        (SubmissionStatus.COMPLETED.name, SubmissionStatus.IN_PROGRESS.name),
        (SubmissionStatus.PENDING.name, SubmissionStatus.APPROVED.name),
        (SubmissionStatus.QUOTE_PROVIDED.name, SubmissionStatus.PENDING.name)
    ]
    
    # For each transition, call submission_service.validate_status_transition
    for current_status, new_status in invalid_transitions:
        # Assert that the result is False for each invalid transition
        with pytest.raises(ValidationException):
            submission_service.validate_status_transition(current_status, new_status)