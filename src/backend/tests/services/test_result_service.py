import pytest  # pytest version: latest
from pytest import fixture, mark, raises  # pytest version: latest
from unittest import mock  # unittest.mock version: standard library
from unittest.mock import MagicMock, patch, Mock  # unittest.mock version: standard library
import io  # io version: standard library
from io import BytesIO  # io version: standard library
import uuid  # uuid version: standard library
from uuid import UUID  # uuid version: standard library
import datetime  # datetime version: standard library

from ..conftest import db_session, test_pharma_user, test_cro_user, mock_file_storage
from app.services.result_service import ResultService
from app.schemas.result import ResultCreate, ResultUpdate, ResultFilter, ResultApproval, ResultDataCreate
from app.constants import ResultStatus, SubmissionStatus
from app.exceptions import ValidationException, ResourceNotFoundException, ResultException

@fixture
def create_test_submission(db_session, pharma_user_id, cro_user_id):
    """Fixture to create a test submission for result testing"""
    from app.models.experiment import Experiment
    from app.models.submission import Submission

    # Create a test experiment in the database
    experiment = Experiment(
        name="Test Experiment",
        type_id=uuid.uuid4(),
        status=SubmissionStatus.IN_PROGRESS.name,
        created_by=pharma_user_id,
        created_at=datetime.datetime.utcnow(),
        updated_at=datetime.datetime.utcnow()
    )
    db_session.add(experiment)
    db_session.flush()

    # Create a test submission linked to the experiment
    submission = Submission(
        experiment_id=experiment.id,
        cro_id=cro_user_id,
        status=SubmissionStatus.IN_PROGRESS.name,
        submitted_at=datetime.datetime.utcnow(),
        updated_at=datetime.datetime.utcnow()
    )
    db_session.add(submission)
    db_session.commit()
    db_session.refresh(submission)

    # Set submission status to IN_PROGRESS
    submission.status = SubmissionStatus.IN_PROGRESS.name
    db_session.commit()

    # Return submission data with ID
    return {"id": submission.id, "experiment_id": experiment.id, "cro_id": cro_user_id}

@fixture
def create_test_result(db_session, test_submission):
    """Fixture to create a test result for testing"""
    from app.models.result import Result

    # Create a test result linked to the test submission
    result = Result(
        submission_id=test_submission["id"],
        status=ResultStatus.UPLOADED.name,
        uploaded_at=datetime.datetime.utcnow()
    )
    db_session.add(result)
    db_session.commit()
    db_session.refresh(result)

    # Set result status to UPLOADED
    result.status = ResultStatus.UPLOADED.name
    db_session.commit()

    # Return result data with ID
    return {"id": result.id, "submission_id": test_submission["id"], "status": ResultStatus.UPLOADED.name}

@fixture
def create_test_molecule(db_session, pharma_user_id):
    """Fixture to create a test molecule for result data testing"""
    from app.models.molecule import Molecule

    # Create a test molecule with SMILES and properties
    molecule = Molecule(
        smiles="CCO",
        created_by=pharma_user_id,
        created_at=datetime.datetime.utcnow()
    )
    db_session.add(molecule)
    db_session.commit()
    db_session.refresh(molecule)

    # Return molecule data with ID
    return {"id": molecule.id, "smiles": "CCO"}

def test_create_result(db_session, test_submission, test_cro_user, mock_file_storage):
    """Test creating a new result"""
    from app.services import submission_service
    from app.worker.tasks.notification_tasks import send_results_notification

    # Mock submission_service.get_submission_by_id to return test submission
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Mock submission_service.update_submission_status
    submission_service.update_submission_status = MagicMock()

    # Mock send_results_uploaded notification function
    send_results_notification = MagicMock()

    # Create ResultCreate object with test submission ID
    result_create = ResultCreate(submission_id=test_submission["id"])

    # Call ResultService.create_result with the data and CRO user ID
    result_service = ResultService()
    result = result_service.create_result(result_create, test_cro_user.id)

    # Assert result is created with correct submission ID
    assert result["submission_id"] == test_submission["id"]

    # Assert result status is set to PENDING
    assert result["status"] == ResultStatus.PENDING.name

    # Assert submission_service.update_submission_status was called
    submission_service.update_submission_status.assert_called_once()

    # Assert send_results_uploaded was called
    send_results_notification.assert_called_once()

def test_create_result_with_files(db_session, test_submission, test_cro_user, mock_file_storage):
    """Test creating a new result with files"""
    from app.services import submission_service
    from app.worker.tasks.notification_tasks import send_results_notification

    # Mock submission_service.get_submission_by_id to return test submission
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Mock submission_service.update_submission_status
    submission_service.update_submission_status = MagicMock()

    # Mock send_results_uploaded notification function
    send_results_notification = MagicMock()

    # Mock upload_result_file to return a file path
    mock_file_storage.upload_result_file = MagicMock(return_value="test_file_path")

    # Create ResultCreate object with test submission ID and files
    result_create = ResultCreate(
        submission_id=test_submission["id"],
        files=[{"file_name": "test_file.txt", "file_path": "test_file_path", "file_size": 1024, "file_type": "text/plain"}]
    )

    # Call ResultService.create_result with the data and CRO user ID
    result_service = ResultService()
    result = result_service.create_result(result_create, test_cro_user.id)

    # Assert result is created with correct submission ID
    assert result["submission_id"] == test_submission["id"]

    # Assert result status is set to UPLOADED (not PENDING because files are included)
    assert result["status"] == ResultStatus.UPLOADED.name

    # Assert upload_result_file was called for each file
    assert mock_file_storage.upload_result_file.call_count == 1

    # Assert submission_service.update_submission_status was called
    submission_service.update_submission_status.assert_called_once()

    # Assert send_results_uploaded was called
    send_results_notification.assert_called_once()

def test_create_result_invalid_submission(db_session, test_cro_user):
    """Test creating a result with invalid submission ID"""
    from app.services import submission_service

    # Mock submission_service.get_submission_by_id to return None
    submission_service.get_submission_by_id = MagicMock(return_value=None)

    # Create ResultCreate object with non-existent submission ID
    result_create = ResultCreate(submission_id=uuid.uuid4())

    # Call ResultService.create_result and expect ResourceNotFoundException
    result_service = ResultService()
    with raises(ResourceNotFoundException) as exc_info:
        result_service.create_result(result_create, test_cro_user.id)

    # Assert error message mentions submission not found
    assert "Submission not found" in str(exc_info.value)

def test_create_result_wrong_cro(db_session, test_submission, test_cro_user):
    """Test creating a result by a CRO not assigned to the submission"""
    from app.services import submission_service

    # Modify test submission to have a different CRO ID
    test_submission["cro_id"] = uuid.uuid4()

    # Mock submission_service.get_submission_by_id to return modified submission
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Create ResultCreate object with test submission ID
    result_create = ResultCreate(submission_id=test_submission["id"])

    # Call ResultService.create_result and expect ValidationException
    result_service = ResultService()
    with raises(ValidationException) as exc_info:
        result_service.create_result(result_create, test_cro_user.id)

    # Assert error message mentions unauthorized CRO
    assert "Unauthorized CRO" in str(exc_info.value)

def test_get_result_by_id(db_session, test_result):
    """Test retrieving a result by ID"""
    # Call ResultService.get_result_by_id with test result ID
    result_service = ResultService()
    result = result_service.get_result_by_id(test_result["id"])

    # Assert returned result has correct ID
    assert result["id"] == test_result["id"]

    # Assert returned result has correct submission ID
    assert result["submission_id"] == test_result["submission_id"]

    # Assert returned result has correct status
    assert result["status"] == test_result["status"]

def test_get_result_by_id_with_data(db_session, test_result, test_molecule):
    """Test retrieving a result by ID with data included"""
    from app.models.result_data import ResultData

    # Add test data point to the result
    data_point = ResultData(
        result_id=test_result["id"],
        molecule_id=test_molecule["id"],
        data_name="BindingAffinity",
        data_value=12.34,
        data_unit="nM"
    )
    db_session.add(data_point)
    db_session.commit()

    # Call ResultService.get_result_by_id with test result ID and include_data=True
    result_service = ResultService()
    result = result_service.get_result_by_id(test_result["id"], include_data=True)

    # Assert returned result has correct ID
    assert result["id"] == test_result["id"]

    # Assert returned result has data_points list
    assert "data_points" in result

    # Assert data point has correct molecule ID and values
    assert len(result["data_points"]) == 1
    assert result["data_points"][0]["molecule_id"] == test_molecule["id"]
    assert result["data_points"][0]["data_name"] == "BindingAffinity"
    assert result["data_points"][0]["data_value"] == 12.34
    assert result["data_points"][0]["data_unit"] == "nM"

def test_update_result(db_session, test_result, test_cro_user):
    """Test updating a result"""
    from app.services import submission_service

    # Mock submission_service.get_submission_by_id to return test submission
    submission_service.get_submission_by_id = MagicMock(return_value={"id": test_result["submission_id"]})

    # Create ResultUpdate object with updated notes
    result_update = ResultUpdate(notes="Updated notes")

    # Call ResultService.update_result with result ID, update data, and CRO user ID
    result_service = ResultService()
    result = result_service.update_result(test_result["id"], result_update, test_cro_user.id)

    # Assert result is updated with new notes
    assert result["notes"] == "Updated notes"

    # Assert other fields remain unchanged
    assert result["submission_id"] == test_result["submission_id"]
    assert result["status"] == test_result["status"]

def test_update_result_not_found(db_session, test_cro_user):
    """Test updating a non-existent result"""
    # Create ResultUpdate object with updated notes
    result_update = ResultUpdate(notes="Updated notes")

    # Call ResultService.update_result with non-existent result ID
    result_service = ResultService()
    with raises(ResourceNotFoundException) as exc_info:
        result_service.update_result(uuid.uuid4(), result_update, test_cro_user.id)

    # Expect ResourceNotFoundException
    assert "Result not found" in str(exc_info.value)

def test_approve_result(db_session, test_result, test_submission, test_pharma_user):
    """Test approving a result"""
    from app.services import submission_service

    # Mock submission_service.get_submission_by_id to return test submission
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Mock submission_service.update_submission_status
    submission_service.update_submission_status = MagicMock()

    # Create ResultApproval object with approved=True
    result_approval = ResultApproval(approved=True)

    # Call ResultService.approve_reject_result with result ID, approval data, and pharma user ID
    result_service = ResultService()
    result = result_service.approve_reject_result(test_result["id"], result_approval, test_pharma_user.id)

    # Assert result status is updated to APPROVED
    assert result["status"] == ResultStatus.APPROVED.name

    # Assert submission_service.update_submission_status was called with COMPLETED status
    submission_service.update_submission_status.assert_called_with(test_submission["id"], SubmissionStatus.COMPLETED.name, test_pharma_user.id)

def test_reject_result(db_session, test_result, test_submission, test_pharma_user):
    """Test rejecting a result"""
    from app.services import submission_service

    # Mock submission_service.get_submission_by_id to return test submission
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Mock submission_service.update_submission_status
    submission_service.update_submission_status = MagicMock()

    # Create ResultApproval object with approved=False and rejection notes
    result_approval = ResultApproval(approved=False, notes="Rejection notes")

    # Call ResultService.approve_reject_result with result ID, approval data, and pharma user ID
    result_service = ResultService()
    result = result_service.approve_reject_result(test_result["id"], result_approval, test_pharma_user.id)

    # Assert result status is updated to REJECTED
    assert result["status"] == ResultStatus.REJECTED.name

    # Assert result notes are updated with rejection reason
    assert result["notes"] == "Rejection notes"

    # Assert submission_service.update_submission_status was called with IN_PROGRESS status
    submission_service.update_submission_status.assert_called_with(test_submission["id"], SubmissionStatus.IN_PROGRESS.name, test_pharma_user.id)

def test_approve_result_unauthorized(db_session, test_result, test_submission, test_cro_user):
    """Test approving a result by unauthorized user"""
    from app.services import submission_service

    # Mock submission_service.get_submission_by_id to return test submission with different creator
    test_submission["created_by"] = uuid.uuid4()
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Create ResultApproval object with approved=True
    result_approval = ResultApproval(approved=True)

    # Call ResultService.approve_reject_result with result ID, approval data, and CRO user ID
    result_service = ResultService()
    with raises(ValidationException) as exc_info:
        result_service.approve_reject_result(test_result["id"], result_approval, test_cro_user.id)

    # Expect ValidationException
    assert "Unauthorized user" in str(exc_info.value)

def test_add_result_file(db_session, test_result, test_submission, test_cro_user, mock_file_storage):
    """Test adding a file to a result"""
    from app.services import submission_service

    # Mock submission_service.get_submission_by_id to return test submission
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Mock upload_result_file to return a file path
    mock_file_storage.upload_result_file = MagicMock(return_value="test_file_path")

    # Create a BytesIO object as test file
    test_file = BytesIO(b"Test file content")

    # Call ResultService.add_result_file with result ID, file object, filename, size, content type, and CRO user ID
    result_service = ResultService()
    result_file = result_service.add_result_file(test_result["id"], test_file, "test_file.txt", 1024, "text/plain", test_cro_user.id)

    # Assert upload_result_file was called with correct parameters
    mock_file_storage.upload_result_file.assert_called_with(test_file, "test_file.txt")

    # Assert file was added to the result
    assert result_file.file_path == "test_file_path"

    # Assert result status is updated to UPLOADED if it was PENDING
    updated_result = db_session.get(Result, test_result["id"])
    assert updated_result.status == ResultStatus.UPLOADED.name

def test_add_result_file_unauthorized(db_session, test_result, test_submission, test_pharma_user, mock_file_storage):
    """Test adding a file to a result by unauthorized user"""
    from app.services import submission_service

    # Mock submission_service.get_submission_by_id to return test submission
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Create a BytesIO object as test file
    test_file = BytesIO(b"Test file content")

    # Call ResultService.add_result_file with result ID, file object, filename, size, content type, and pharma user ID
    result_service = ResultService()
    with raises(ValidationException) as exc_info:
        result_service.add_result_file(test_result["id"], test_file, "test_file.txt", 1024, "text/plain", test_pharma_user.id)

    # Expect ValidationException
    assert "Unauthorized user" in str(exc_info.value)

def test_get_result_file(db_session, test_result, test_submission, test_pharma_user, mock_file_storage):
    """Test retrieving a result file"""
    from app.services import submission_service
    from app.models.result_file import ResultFile

    # Add a test file to the result
    result_file = ResultFile(result_id=test_result["id"], file_name="test_file.txt", file_path="test_file_path", file_size=1024, file_type="text/plain")
    db_session.add(result_file)
    db_session.commit()

    # Mock submission_service.get_submission_by_id to return test submission
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Mock get_result_file to return file data, filename, and content type
    mock_file_storage.get_result_file = MagicMock(return_value=(b"Test file content", "test_file.txt", "text/plain"))

    # Call ResultService.get_result_file with file ID and pharma user ID
    result_service = ResultService()
    file_data, filename, content_type = result_service.get_result_file(result_file.id, test_pharma_user.id)

    # Assert get_result_file was called with correct parameters
    mock_file_storage.get_result_file.assert_called_with("test_file_path")

    # Assert returned tuple contains file data, filename, and content type
    assert file_data == b"Test file content"
    assert filename == "test_file.txt"
    assert content_type == "text/plain"

def test_get_result_file_url(db_session, test_result, test_submission, test_pharma_user, mock_file_storage):
    """Test generating a URL for a result file"""
    from app.services import submission_service
    from app.models.result_file import ResultFile

    # Add a test file to the result
    result_file = ResultFile(result_id=test_result["id"], file_name="test_file.txt", file_path="test_file_path", file_size=1024, file_type="text/plain")
    db_session.add(result_file)
    db_session.commit()

    # Mock submission_service.get_submission_by_id to return test submission
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Mock get_result_file_url to return a presigned URL
    mock_file_storage.get_result_file_url = MagicMock(return_value="http://test_url")

    # Call ResultService.get_result_file_url with file ID, pharma user ID, and expiration
    result_service = ResultService()
    url = result_service.get_result_file_url(result_file.id, test_pharma_user.id, 3600)

    # Assert get_result_file_url was called with correct parameters
    mock_file_storage.get_result_file_url.assert_called_with("test_file_path", expires=3600)

    # Assert returned URL matches expected URL
    assert url == "http://test_url"

def test_delete_result_file(db_session, test_result, test_submission, test_cro_user, mock_file_storage):
    """Test deleting a result file"""
    from app.services import submission_service
    from app.models.result_file import ResultFile

    # Add a test file to the result
    result_file = ResultFile(result_id=test_result["id"], file_name="test_file.txt", file_path="test_file_path", file_size=1024, file_type="text/plain")
    db_session.add(result_file)
    db_session.commit()

    # Mock submission_service.get_submission_by_id to return test submission
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Mock delete_file to return True
    mock_file_storage.delete_file = MagicMock(return_value=True)

    # Call ResultService.delete_result_file with file ID and CRO user ID
    result_service = ResultService()
    result = result_service.delete_result_file(result_file.id, test_cro_user.id)

    # Assert delete_file was called with correct parameters
    mock_file_storage.delete_file.assert_called_with("test_file_path")

    # Assert file was removed from the database
    assert db_session.get(ResultFile, result_file.id) is None

def test_delete_result_file_unauthorized(db_session, test_result, test_submission, test_pharma_user, mock_file_storage):
    """Test deleting a result file by unauthorized user"""
    from app.services import submission_service
    from app.models.result_file import ResultFile

    # Add a test file to the result
    result_file = ResultFile(result_id=test_result["id"], file_name="test_file.txt", file_path="test_file_path", file_size=1024, file_type="text/plain")
    db_session.add(result_file)
    db_session.commit()

    # Mock submission_service.get_submission_by_id to return test submission
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Call ResultService.delete_result_file with file ID and pharma user ID
    result_service = ResultService()
    with raises(ValidationException) as exc_info:
        result_service.delete_result_file(result_file.id, test_pharma_user.id)

    # Expect ValidationException
    assert "Unauthorized user" in str(exc_info.value)

def test_add_result_data(db_session, test_result, test_submission, test_molecule, test_cro_user):
    """Test adding structured data to a result"""
    from app.services import submission_service

    # Mock submission_service.get_submission_by_id to return test submission
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Create ResultDataCreate object with molecule ID, data name, value, and unit
    result_data_create = ResultDataCreate(
        molecule_id=test_molecule["id"],
        data_name="BindingAffinity",
        data_value=12.34,
        data_unit="nM"
    )

    # Call ResultService.add_result_data with result ID, data item, and CRO user ID
    result_service = ResultService()
    result_data = result_service.add_result_data(test_result["id"], result_data_create, test_cro_user.id)

    # Assert data was added to the result
    assert result_data.molecule_id == test_molecule["id"]
    assert result_data.data_name == "BindingAffinity"
    assert result_data.data_value == 12.34
    assert result_data.data_unit == "nM"

    # Assert result status is updated to UPLOADED if it was PENDING
    updated_result = db_session.get(Result, test_result["id"])
    assert updated_result.status == ResultStatus.UPLOADED.name

def test_get_result_data(db_session, test_result, test_submission, test_molecule, test_pharma_user):
    """Test retrieving structured data for a result"""
    from app.services import submission_service
    from app.models.result_data import ResultData

    # Add test data point to the result
    data_point = ResultData(
        result_id=test_result["id"],
        molecule_id=test_molecule["id"],
        data_name="BindingAffinity",
        data_value=12.34,
        data_unit="nM"
    )
    db_session.add(data_point)
    db_session.commit()

    # Mock submission_service.get_submission_by_id to return test submission
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Call ResultService.get_result_data with result ID and pharma user ID
    result_service = ResultService()
    data = result_service.get_result_data(test_result["id"], test_pharma_user.id)

    # Assert returned data list contains the test data point
    assert len(data) == 1
    assert data[0].molecule_id == test_molecule["id"]
    assert data[0].data_name == "BindingAffinity"
    assert data[0].data_value == 12.34
    assert data[0].data_unit == "nM"

def test_get_results_by_submission(db_session, test_result, test_submission, test_pharma_user):
    """Test retrieving results for a submission"""
    from app.services import submission_service

    # Mock submission_service.get_submission_by_id to return test submission
    submission_service.get_submission_by_id = MagicMock(return_value=test_submission)

    # Call ResultService.get_results_by_submission with submission ID, pharma user ID, skip, and limit
    result_service = ResultService()
    results = result_service.get_results_by_submission(test_result["submission_id"], test_pharma_user.id, skip=0, limit=10)

    # Assert returned dictionary contains items list with the test result
    assert "items" in results
    assert len(results["items"]) == 1
    assert results["items"][0]["id"] == test_result["id"]

    # Assert total count is 1
    assert results["total"] == 1

    # Assert page and size match expected values
    assert results["page"] == 1
    assert results["size"] == 10

def test_get_results(db_session, test_result, test_pharma_user):
    """Test retrieving results with filtering"""
    # Create ResultFilter object with status filter
    result_filter = ResultFilter(status=ResultStatus.UPLOADED.name)

    # Call ResultService.get_results with filter, pharma user ID, skip, and limit
    result_service = ResultService()
    results = result_service.get_results(result_filter, test_pharma_user.id, skip=0, limit=10)

    # Assert returned dictionary contains items list with the test result
    assert "items" in results
    assert len(results["items"]) == 1
    assert results["items"][0]["id"] == test_result["id"]

    # Assert total count is 1
    assert results["total"] == 1

    # Assert page and size match expected values
    assert results["page"] == 1
    assert results["size"] == 10

def test_validate_status_transition():
    """Test validating result status transitions"""
    # Create ResultService instance
    result_service = ResultService()

    # Test valid transitions: PENDING to UPLOADED, UPLOADED to APPROVED, etc.
    assert result_service.validate_status_transition(ResultStatus.PENDING.name, ResultStatus.UPLOADED.name) is None
    assert result_service.validate_status_transition(ResultStatus.UPLOADED.name, ResultStatus.APPROVED.name) is None
    assert result_service.validate_status_transition(ResultStatus.UPLOADED.name, ResultStatus.REJECTED.name) is None

    # Test invalid transitions: APPROVED to UPLOADED, REJECTED to PENDING, etc.
    with raises(ValidationException):
        result_service.validate_status_transition(ResultStatus.APPROVED.name, ResultStatus.UPLOADED.name)
    with raises(ValidationException):
        result_service.validate_status_transition(ResultStatus.REJECTED.name, ResultStatus.PENDING.name)