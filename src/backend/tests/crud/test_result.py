import pytest
from datetime import datetime, timedelta
import uuid
import random

from app.crud.crud_result import result
from app.models.result import Result, ResultStatus
from app.models.result_file import ResultFile
from app.models.result_data import ResultData
from app.schemas.result import ResultCreate, ResultUpdate


def test_create_result(db_session, test_cro_user):
    """Test creating a new result"""
    # Create a test submission for the result
    experiment_id = uuid.uuid4()
    submission = create_test_submission(db_session, test_cro_user.id, experiment_id)
    
    # Create result
    result_in = ResultCreate(
        submission_id=submission.id,
        status=ResultStatus.PENDING.name,
        notes="Test result notes"
    )
    
    # Call create function
    result_obj = result.create(db_session, result_in)
    
    # Verify result was created with correct data
    assert result_obj.submission_id == submission.id
    assert result_obj.status == ResultStatus.PENDING.name
    assert result_obj.notes == "Test result notes"
    assert result_obj.uploaded_at is not None
    assert (datetime.utcnow() - result_obj.uploaded_at).total_seconds() < 10  # Created within last 10 seconds
    assert result_obj.approved_at is None


def test_get_result(db_session, test_cro_user):
    """Test retrieving a result by ID"""
    # Create a test submission and result
    experiment_id = uuid.uuid4()
    submission = create_test_submission(db_session, test_cro_user.id, experiment_id)
    result_obj = create_test_result(db_session, submission.id)
    
    # Get the result by ID
    retrieved_result = result.get(db_session, result_obj.id)
    
    # Verify result was retrieved correctly
    assert retrieved_result is not None
    assert retrieved_result.id == result_obj.id
    assert retrieved_result.submission_id == submission.id
    
    # Test with non-existent ID
    non_existent_id = uuid.uuid4()
    assert result.get(db_session, non_existent_id) is None


def test_update_result(db_session, test_cro_user):
    """Test updating a result"""
    # Create a test submission and result
    experiment_id = uuid.uuid4()
    submission = create_test_submission(db_session, test_cro_user.id, experiment_id)
    result_obj = create_test_result(db_session, submission.id)
    
    # Get the result object to update
    db_obj = result.get(db_session, result_obj.id)
    
    # Update the result
    result_update = ResultUpdate(
        status=ResultStatus.UPLOADED.name,
        notes="Updated test notes"
    )
    
    updated_result = result.update(db_session, db_obj, result_update)
    
    # Verify result was updated correctly
    assert updated_result.status == ResultStatus.UPLOADED.name
    assert updated_result.notes == "Updated test notes"
    
    # Test updating with dictionary instead of schema
    db_obj = result.get(db_session, result_obj.id)
    dict_update = {
        "status": ResultStatus.APPROVED.name,
        "notes": "Dict update notes"
    }
    
    updated_result = result.update(db_session, db_obj, dict_update)
    
    # Verify result was updated correctly
    assert updated_result.status == ResultStatus.APPROVED.name
    assert updated_result.notes == "Dict update notes"


def test_get_with_details(db_session, test_cro_user, create_test_molecule):
    """Test retrieving a result with all related details"""
    # Create a test submission and result
    experiment_id = uuid.uuid4()
    submission = create_test_submission(db_session, test_cro_user.id, experiment_id)
    result_obj = create_test_result(db_session, submission.id)
    
    # Add test files
    file1 = result.add_file(db_session, result_obj.id, "test1.pdf", "/path/to/test1.pdf", 1024, "application/pdf")
    file2 = result.add_file(db_session, result_obj.id, "test2.csv", "/path/to/test2.csv", 2048, "text/csv")
    
    # Add test data points
    molecule1 = create_test_molecule()
    molecule2 = create_test_molecule()
    data1 = result.add_data(db_session, result_obj.id, molecule1.id, "IC50", 12.5, "nM")
    data2 = result.add_data(db_session, result_obj.id, molecule2.id, "Solubility", 0.25, "mg/mL")
    
    # Get result with details
    detailed_result = result.get_with_details(db_session, result_obj.id)
    
    # Verify result and related data
    assert detailed_result is not None
    assert detailed_result.id == result_obj.id
    assert detailed_result.submission is not None
    assert len(detailed_result.files) == 2
    assert len(detailed_result.data_points) == 2
    
    # Verify file data
    file_paths = [f.file_path for f in detailed_result.files]
    assert "/path/to/test1.pdf" in file_paths
    assert "/path/to/test2.csv" in file_paths
    
    # Verify data points
    data_names = [d.data_name for d in detailed_result.data_points]
    assert "IC50" in data_names
    assert "Solubility" in data_names


def test_get_with_files(db_session, test_cro_user):
    """Test retrieving a result with its associated files"""
    # Create a test submission and result
    experiment_id = uuid.uuid4()
    submission = create_test_submission(db_session, test_cro_user.id, experiment_id)
    result_obj = create_test_result(db_session, submission.id)
    
    # Add test files
    file1 = result.add_file(db_session, result_obj.id, "test1.pdf", "/path/to/test1.pdf", 1024, "application/pdf")
    file2 = result.add_file(db_session, result_obj.id, "test2.csv", "/path/to/test2.csv", 2048, "text/csv")
    
    # Get result with files
    result_with_files = result.get_with_files(db_session, result_obj.id)
    
    # Verify result and files
    assert result_with_files is not None
    assert len(result_with_files.files) == 2
    
    file_names = [f.file_name for f in result_with_files.files]
    file_paths = [f.file_path for f in result_with_files.files]
    file_sizes = [f.file_size for f in result_with_files.files]
    file_types = [f.file_type for f in result_with_files.files]
    
    assert "test1.pdf" in file_names
    assert "test2.csv" in file_names
    assert "/path/to/test1.pdf" in file_paths
    assert "/path/to/test2.csv" in file_paths
    assert 1024 in file_sizes
    assert 2048 in file_sizes
    assert "application/pdf" in file_types
    assert "text/csv" in file_types


def test_get_with_data(db_session, test_cro_user, create_test_molecule):
    """Test retrieving a result with its associated data points"""
    # Create a test submission and result
    experiment_id = uuid.uuid4()
    submission = create_test_submission(db_session, test_cro_user.id, experiment_id)
    result_obj = create_test_result(db_session, submission.id)
    
    # Create test molecules
    molecule1 = create_test_molecule()
    molecule2 = create_test_molecule()
    
    # Add test data points
    data1 = result.add_data(db_session, result_obj.id, molecule1.id, "IC50", 12.5, "nM")
    data2 = result.add_data(db_session, result_obj.id, molecule2.id, "Solubility", 0.25, "mg/mL")
    
    # Get result with data
    result_with_data = result.get_with_data(db_session, result_obj.id)
    
    # Verify result and data points
    assert result_with_data is not None
    assert len(result_with_data.data_points) == 2
    
    data_names = [d.data_name for d in result_with_data.data_points]
    data_values = [d.data_value for d in result_with_data.data_points]
    data_units = [d.data_unit for d in result_with_data.data_points]
    
    assert "IC50" in data_names
    assert "Solubility" in data_names
    assert 12.5 in data_values
    assert 0.25 in data_values
    assert "nM" in data_units
    assert "mg/mL" in data_units


def test_get_by_submission(db_session, test_cro_user):
    """Test retrieving results for a specific submission"""
    # Create a test submission
    experiment_id = uuid.uuid4()
    submission = create_test_submission(db_session, test_cro_user.id, experiment_id)
    
    # Create multiple results for the submission
    result1 = create_test_result(db_session, submission.id, ResultStatus.PENDING.name)
    result2 = create_test_result(db_session, submission.id, ResultStatus.UPLOADED.name)
    result3 = create_test_result(db_session, submission.id, ResultStatus.APPROVED.name)
    
    # Get results for the submission
    submission_results = result.get_by_submission(db_session, submission.id)
    
    # Verify results
    assert len(submission_results) == 3
    result_ids = [r.id for r in submission_results]
    assert result1.id in result_ids
    assert result2.id in result_ids
    assert result3.id in result_ids
    
    # Test pagination
    paginated_results = result.get_by_submission(db_session, submission.id, skip=1, limit=1)
    assert len(paginated_results) == 1
    
    # Test with non-existent submission ID
    non_existent_id = uuid.uuid4()
    assert len(result.get_by_submission(db_session, non_existent_id)) == 0


def test_get_by_status(db_session, test_cro_user):
    """Test retrieving results with a specific status"""
    # Create test submissions
    experiment_id1 = uuid.uuid4()
    experiment_id2 = uuid.uuid4()
    submission1 = create_test_submission(db_session, test_cro_user.id, experiment_id1)
    submission2 = create_test_submission(db_session, test_cro_user.id, experiment_id2)
    
    # Create results with different statuses
    result1 = create_test_result(db_session, submission1.id, ResultStatus.PENDING.name)
    result2 = create_test_result(db_session, submission1.id, ResultStatus.UPLOADED.name)
    result3 = create_test_result(db_session, submission2.id, ResultStatus.UPLOADED.name)
    result4 = create_test_result(db_session, submission2.id, ResultStatus.APPROVED.name)
    
    # Get results by status
    pending_results = result.get_by_status(db_session, ResultStatus.PENDING.name)
    uploaded_results = result.get_by_status(db_session, ResultStatus.UPLOADED.name)
    approved_results = result.get_by_status(db_session, ResultStatus.APPROVED.name)
    
    # Verify results
    assert len(pending_results) == 1
    assert pending_results[0].id == result1.id
    
    assert len(uploaded_results) == 2
    uploaded_ids = [r.id for r in uploaded_results]
    assert result2.id in uploaded_ids
    assert result3.id in uploaded_ids
    
    assert len(approved_results) == 1
    assert approved_results[0].id == result4.id
    
    # Test pagination
    paginated_results = result.get_by_status(db_session, ResultStatus.UPLOADED.name, skip=1, limit=1)
    assert len(paginated_results) == 1


def test_update_status(db_session, test_cro_user):
    """Test updating the status of a result"""
    # Create a test submission and result with initial status PENDING
    experiment_id = uuid.uuid4()
    submission = create_test_submission(db_session, test_cro_user.id, experiment_id)
    result_obj = create_test_result(db_session, submission.id, ResultStatus.PENDING.name)
    
    # Update status to UPLOADED
    updated_result = result.update_status(db_session, result_obj.id, ResultStatus.UPLOADED.name)
    
    # Verify status was updated
    assert updated_result.status == ResultStatus.UPLOADED.name
    assert updated_result.approved_at is None
    
    # Update status to APPROVED (should set approved_at)
    updated_result = result.update_status(db_session, result_obj.id, ResultStatus.APPROVED.name)
    
    # Verify status and approved_at
    assert updated_result.status == ResultStatus.APPROVED.name
    assert updated_result.approved_at is not None
    assert (datetime.utcnow() - updated_result.approved_at).total_seconds() < 10  # Set within last 10 seconds
    
    # Test with non-existent result ID
    non_existent_id = uuid.uuid4()
    assert result.update_status(db_session, non_existent_id, ResultStatus.UPLOADED.name) is None


def test_add_file(db_session, test_cro_user):
    """Test adding a file to a result"""
    # Create a test submission and result
    experiment_id = uuid.uuid4()
    submission = create_test_submission(db_session, test_cro_user.id, experiment_id)
    result_obj = create_test_result(db_session, submission.id)
    
    # Add a file
    file = result.add_file(
        db_session, 
        result_obj.id, 
        "test.pdf", 
        "/path/to/test.pdf", 
        1024, 
        "application/pdf"
    )
    
    # Verify file was added correctly
    assert file.result_id == result_obj.id
    assert file.file_name == "test.pdf"
    assert file.file_path == "/path/to/test.pdf"
    assert file.file_size == 1024
    assert file.file_type == "application/pdf"
    assert file.uploaded_at is not None
    assert (datetime.utcnow() - file.uploaded_at).total_seconds() < 10  # Created within last 10 seconds
    
    # Add multiple files and verify they're associated with the result
    file2 = result.add_file(db_session, result_obj.id, "test2.xlsx", "/path/to/test2.xlsx", 2048, "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    file3 = result.add_file(db_session, result_obj.id, "test3.csv", "/path/to/test3.csv", 512, "text/csv")
    
    # Get result with files
    result_with_files = result.get_with_files(db_session, result_obj.id)
    
    # Verify all files are associated with the result
    assert len(result_with_files.files) == 3


def test_get_file(db_session, test_cro_user):
    """Test retrieving a result file by ID"""
    # Create a test submission and result
    experiment_id = uuid.uuid4()
    submission = create_test_submission(db_session, test_cro_user.id, experiment_id)
    result_obj = create_test_result(db_session, submission.id)
    
    # Add a file
    file = result.add_file(
        db_session, 
        result_obj.id, 
        "test.pdf", 
        "/path/to/test.pdf", 
        1024, 
        "application/pdf"
    )
    
    # Get the file
    retrieved_file = result.get_file(db_session, file.id)
    
    # Verify file was retrieved correctly
    assert retrieved_file is not None
    assert retrieved_file.id == file.id
    assert retrieved_file.result_id == result_obj.id
    assert retrieved_file.file_name == "test.pdf"
    assert retrieved_file.file_path == "/path/to/test.pdf"
    assert retrieved_file.file_size == 1024
    assert retrieved_file.file_type == "application/pdf"
    
    # Test with non-existent ID
    non_existent_id = uuid.uuid4()
    assert result.get_file(db_session, non_existent_id) is None


def test_delete_file(db_session, test_cro_user):
    """Test deleting a result file"""
    # Create a test submission and result
    experiment_id = uuid.uuid4()
    submission = create_test_submission(db_session, test_cro_user.id, experiment_id)
    result_obj = create_test_result(db_session, submission.id)
    
    # Add a file
    file = result.add_file(
        db_session, 
        result_obj.id, 
        "test.pdf", 
        "/path/to/test.pdf", 
        1024, 
        "application/pdf"
    )
    
    # Delete the file
    delete_success = result.delete_file(db_session, file.id)
    
    # Verify deletion was successful
    assert delete_success is True
    
    # Verify file no longer exists
    assert result.get_file(db_session, file.id) is None
    
    # Test with non-existent ID
    non_existent_id = uuid.uuid4()
    assert result.delete_file(db_session, non_existent_id) is False


def test_add_data(db_session, test_cro_user, create_test_molecule):
    """Test adding structured data to a result for a specific molecule"""
    # Create a test submission and result
    experiment_id = uuid.uuid4()
    submission = create_test_submission(db_session, test_cro_user.id, experiment_id)
    result_obj = create_test_result(db_session, submission.id)
    
    # Create a test molecule
    molecule = create_test_molecule()
    
    # Add data
    data = result.add_data(
        db_session,
        result_obj.id,
        molecule.id,
        "IC50",
        12.5,
        "nM"
    )
    
    # Verify data was added correctly
    assert data.result_id == result_obj.id
    assert data.molecule_id == molecule.id
    assert data.data_name == "IC50"
    assert data.data_value == 12.5
    assert data.data_unit == "nM"
    
    # Add multiple data points and verify they're associated with the result
    data2 = result.add_data(db_session, result_obj.id, molecule.id, "Solubility", 0.25, "mg/mL")
    data3 = result.add_data(db_session, result_obj.id, molecule.id, "LogP", 2.3, None)
    
    # Get result data
    result_data = result.get_data(db_session, result_obj.id)
    
    # Verify all data points are associated with the result
    assert len(result_data) == 3


def test_get_data(db_session, test_cro_user, create_test_molecule):
    """Test retrieving all structured data for a result"""
    # Create a test submission and result
    experiment_id = uuid.uuid4()
    submission = create_test_submission(db_session, test_cro_user.id, experiment_id)
    result_obj = create_test_result(db_session, submission.id)
    
    # Create test molecules
    molecule1 = create_test_molecule()
    molecule2 = create_test_molecule()
    
    # Add test data points
    data1 = result.add_data(db_session, result_obj.id, molecule1.id, "IC50", 12.5, "nM")
    data2 = result.add_data(db_session, result_obj.id, molecule1.id, "LogP", 2.3, None)
    data3 = result.add_data(db_session, result_obj.id, molecule2.id, "Solubility", 0.25, "mg/mL")
    
    # Get all result data
    result_data = result.get_data(db_session, result_obj.id)
    
    # Verify data was retrieved correctly
    assert len(result_data) == 3
    
    data_names = [d.data_name for d in result_data]
    data_values = [d.data_value for d in result_data]
    data_units = [d.data_unit for d in result_data]
    molecule_ids = [d.molecule_id for d in result_data]
    
    assert "IC50" in data_names
    assert "LogP" in data_names
    assert "Solubility" in data_names
    assert 12.5 in data_values
    assert 2.3 in data_values
    assert 0.25 in data_values
    assert "nM" in data_units
    assert None in data_units
    assert "mg/mL" in data_units
    assert molecule1.id in molecule_ids
    assert molecule2.id in molecule_ids
    
    # Test with non-existent result ID
    non_existent_id = uuid.uuid4()
    assert len(result.get_data(db_session, non_existent_id)) == 0


def test_get_data_by_molecule(db_session, test_cro_user, create_test_molecule):
    """Test retrieving structured data for a result filtered by molecule"""
    # Create a test submission and result
    experiment_id = uuid.uuid4()
    submission = create_test_submission(db_session, test_cro_user.id, experiment_id)
    result_obj = create_test_result(db_session, submission.id)
    
    # Create test molecules
    molecule1 = create_test_molecule()
    molecule2 = create_test_molecule()
    molecule3 = create_test_molecule()
    
    # Add test data points
    data1 = result.add_data(db_session, result_obj.id, molecule1.id, "IC50", 12.5, "nM")
    data2 = result.add_data(db_session, result_obj.id, molecule1.id, "LogP", 2.3, None)
    data3 = result.add_data(db_session, result_obj.id, molecule2.id, "Solubility", 0.25, "mg/mL")
    data4 = result.add_data(db_session, result_obj.id, molecule2.id, "pKa", 4.5, None)
    data5 = result.add_data(db_session, result_obj.id, molecule3.id, "MW", 250.3, "g/mol")
    
    # Get data for molecule1
    molecule1_data = result.get_data_by_molecule(db_session, result_obj.id, molecule1.id)
    
    # Verify data for molecule1
    assert len(molecule1_data) == 2
    data_names = [d.data_name for d in molecule1_data]
    assert "IC50" in data_names
    assert "LogP" in data_names
    
    # Get data for molecule2
    molecule2_data = result.get_data_by_molecule(db_session, result_obj.id, molecule2.id)
    
    # Verify data for molecule2
    assert len(molecule2_data) == 2
    data_names = [d.data_name for d in molecule2_data]
    assert "Solubility" in data_names
    assert "pKa" in data_names
    
    # Get data for molecule3
    molecule3_data = result.get_data_by_molecule(db_session, result_obj.id, molecule3.id)
    
    # Verify data for molecule3
    assert len(molecule3_data) == 1
    assert molecule3_data[0].data_name == "MW"
    
    # Test with non-existent molecule ID
    non_existent_id = uuid.uuid4()
    assert len(result.get_data_by_molecule(db_session, result_obj.id, non_existent_id)) == 0


def test_filter_results(db_session, test_cro_user, test_pharma_user):
    """Test filtering results based on various criteria"""
    # Create test submissions with different attributes
    experiment_id1 = uuid.uuid4()
    experiment_id2 = uuid.uuid4()
    
    submission1 = create_test_submission(db_session, test_cro_user.id, experiment_id1)
    submission2 = create_test_submission(db_session, test_cro_user.id, experiment_id2)
    
    # Create results with different statuses and dates
    result1 = create_test_result(db_session, submission1.id, ResultStatus.PENDING.name, 
                                datetime.utcnow() - timedelta(days=5))
    result2 = create_test_result(db_session, submission1.id, ResultStatus.UPLOADED.name, 
                                datetime.utcnow() - timedelta(days=3))
    result3 = create_test_result(db_session, submission2.id, ResultStatus.UPLOADED.name, 
                                datetime.utcnow() - timedelta(days=2))
    result4 = create_test_result(db_session, submission2.id, ResultStatus.APPROVED.name, 
                                datetime.utcnow() - timedelta(days=1))
    
    # Test filtering by submission_id
    submission_filter = {"submission_id": submission1.id}
    filtered_results, count = result.filter_results(db_session, submission_filter)
    
    assert count == 2
    result_ids = [r.id for r in filtered_results]
    assert result1.id in result_ids
    assert result2.id in result_ids
    
    # Test filtering by status
    status_filter = {"status": ResultStatus.UPLOADED.name}
    filtered_results, count = result.filter_results(db_session, status_filter)
    
    assert count == 2
    result_ids = [r.id for r in filtered_results]
    assert result2.id in result_ids
    assert result3.id in result_ids
    
    # Test filtering by date range
    date_filter = {
        "uploaded_after": datetime.utcnow() - timedelta(days=4),
        "uploaded_before": datetime.utcnow() - timedelta(days=1, hours=12)
    }
    filtered_results, count = result.filter_results(db_session, date_filter)
    
    assert count == 2
    result_ids = [r.id for r in filtered_results]
    assert result2.id in result_ids
    assert result3.id in result_ids
    
    # Test pagination
    pagination_filter = {
        "skip": 1,
        "limit": 2
    }
    filtered_results, count = result.filter_results(db_session, pagination_filter)
    
    assert len(filtered_results) == 2
    assert count == 4  # Total count is still 4


# Helper functions for test data setup

def create_test_submission(db_session, cro_id, experiment_id, user_id=None):
    """Create a test submission in the database"""
    from app.models.submission import Submission
    
    # Create submission
    submission = Submission(
        experiment_id=experiment_id,
        cro_id=cro_id,
        status="PENDING",
        submitted_at=datetime.utcnow(),
        updated_at=datetime.utcnow()
    )
    
    db_session.add(submission)
    db_session.commit()
    db_session.refresh(submission)
    
    return submission


def create_test_result(db_session, submission_id, status=None, uploaded_at=None):
    """Create a test result in the database"""
    
    if status is None:
        status = ResultStatus.PENDING.name
    
    if uploaded_at is None:
        uploaded_at = datetime.utcnow()
    
    # Create result
    result_obj = Result(
        submission_id=submission_id,
        status=status,
        uploaded_at=uploaded_at,
        notes=f"Test notes for {uuid.uuid4()}"
    )
    
    db_session.add(result_obj)
    db_session.commit()
    db_session.refresh(result_obj)
    
    return result_obj