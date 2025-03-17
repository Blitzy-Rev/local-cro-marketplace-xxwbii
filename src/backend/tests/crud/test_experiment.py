import pytest
import uuid
from typing import List, Dict, Any

from ../../app/crud/crud_experiment import experiment
from ../../app/models/experiment import Experiment
from ../../app/models/experiment_parameter import ExperimentParameter
from ../../app/models/experiment_type import ExperimentType
from ../../app/schemas/experiment import ExperimentCreate, ExperimentUpdate
from ../../app/constants import ExperimentStatus
from ../conftest import db_session, test_pharma_user, create_test_experiment, create_test_molecule


def test_get_by_name(db_session, create_test_experiment, test_pharma_user):
    """Test retrieving an experiment by its name and creator ID"""
    # Create a test experiment with a specific name
    test_name = "Test Experiment for get_by_name"
    test_exp = create_test_experiment(
        name=test_name,
        created_by=test_pharma_user.id,
        experiment_type_id=1,  # Assuming we have this experiment type
        status=ExperimentStatus.DRAFT.name
    )
    
    # Retrieve the experiment using get_by_name function
    retrieved_exp = experiment.get_by_name(db_session, test_name, test_pharma_user.id)
    
    # Assert that the retrieved experiment matches the created one
    assert retrieved_exp.id == test_exp.id
    assert retrieved_exp.name == test_name
    assert retrieved_exp.created_by == test_pharma_user.id
    
    # Test with a non-existent name
    non_existent = experiment.get_by_name(db_session, "Non-existent Experiment", test_pharma_user.id)
    assert non_existent is None


def test_create_experiment(db_session, test_pharma_user):
    """Test creating a new experiment with parameters"""
    # Create an experiment type for testing
    exp_type = ExperimentType(name="Binding Assay", description="Test binding assay", category="Binding")
    db_session.add(exp_type)
    db_session.flush()
    
    # Create experiment data with parameters
    experiment_data = ExperimentCreate(
        name="Test Create Experiment",
        type_id=exp_type.id,
        status=ExperimentStatus.DRAFT.name,
        parameters=[
            {"parameter_name": "concentration", "parameter_value": "10 nM"},
            {"parameter_name": "temperature", "parameter_value": "25°C"}
        ]
    )
    
    # Call the create function with the experiment data
    new_exp = experiment.create(db_session, experiment_data, test_pharma_user.id)
    
    # Assert that the experiment was created with the correct attributes
    assert new_exp.name == "Test Create Experiment"
    assert new_exp.type_id == exp_type.id
    assert new_exp.status == ExperimentStatus.DRAFT
    assert new_exp.created_by == test_pharma_user.id
    
    # Get parameters and check they were created correctly
    params = db_session.query(ExperimentParameter).filter(
        ExperimentParameter.experiment_id == new_exp.id
    ).all()
    
    assert len(params) == 2
    param_dict = {p.parameter_name: p.parameter_value for p in params}
    assert param_dict["concentration"] == "10 nM"
    assert param_dict["temperature"] == "25°C"
    
    # Test creating an experiment without parameters
    experiment_data_no_params = ExperimentCreate(
        name="Test Create No Params",
        type_id=exp_type.id,
        status=ExperimentStatus.DRAFT.name,
        parameters=[]
    )
    
    new_exp_no_params = experiment.create(db_session, experiment_data_no_params, test_pharma_user.id)
    
    # Assert that the experiment is created without parameters
    assert new_exp_no_params.name == "Test Create No Params"
    params_no_params = db_session.query(ExperimentParameter).filter(
        ExperimentParameter.experiment_id == new_exp_no_params.id
    ).all()
    assert len(params_no_params) == 0


def test_update_experiment(db_session, create_test_experiment):
    """Test updating an existing experiment and its parameters"""
    # Create a test experiment with parameters
    exp_type = ExperimentType(name="Update Test", description="Test update", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    test_exp = create_test_experiment(
        name="Test Update Experiment",
        created_by=1,  # Assuming user ID 1 exists
        experiment_type_id=exp_type.id,
        parameters={
            "initial_param": "initial_value"
        }
    )
    
    # Create update data with modified name and parameters
    update_data = ExperimentUpdate(
        name="Updated Experiment Name",
        parameters=[
            {"parameter_name": "initial_param", "parameter_value": "updated_value"},
            {"parameter_name": "new_param", "parameter_value": "new_value"}
        ]
    )
    
    # Call the update function with the experiment and update data
    updated_exp = experiment.update(db_session, test_exp, update_data)
    
    # Assert that the experiment attributes were updated correctly
    assert updated_exp.id == test_exp.id
    assert updated_exp.name == "Updated Experiment Name"
    
    # Check parameters were updated correctly
    params = db_session.query(ExperimentParameter).filter(
        ExperimentParameter.experiment_id == updated_exp.id
    ).all()
    
    assert len(params) == 2
    param_dict = {p.parameter_name: p.parameter_value for p in params}
    assert param_dict["initial_param"] == "updated_value"
    assert param_dict["new_param"] == "new_value"
    
    # Test updating only specific fields
    partial_update_data = ExperimentUpdate(
        status=ExperimentStatus.QUEUED.name
    )
    
    partially_updated_exp = experiment.update(db_session, updated_exp, partial_update_data)
    
    # Assert that only the specified fields were updated
    assert partially_updated_exp.id == test_exp.id
    assert partially_updated_exp.name == "Updated Experiment Name"  # Unchanged
    assert partially_updated_exp.status == ExperimentStatus.QUEUED  # Updated


def test_get_with_parameters(db_session, create_test_experiment):
    """Test retrieving an experiment with its parameters"""
    # Create a test experiment with parameters
    exp_type = ExperimentType(name="Param Test", description="Test params", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    test_exp = create_test_experiment(
        name="Test Get Parameters",
        created_by=1,  # Assuming user ID 1 exists
        experiment_type_id=exp_type.id,
        parameters={
            "param1": "value1",
            "param2": "value2"
        }
    )
    
    # Retrieve the experiment with parameters using get_with_parameters function
    retrieved_exp = experiment.get_with_parameters(db_session, test_exp.id)
    
    # Assert that the experiment has the correct number of parameters
    assert retrieved_exp.id == test_exp.id
    assert hasattr(retrieved_exp, "parameters")
    assert len(retrieved_exp.parameters) == 2
    
    # Assert that the parameter values match the expected values
    param_dict = {p.parameter_name: p.parameter_value for p in retrieved_exp.parameters}
    assert param_dict["param1"] == "value1"
    assert param_dict["param2"] == "value2"
    
    # Test with a non-existent experiment ID
    non_existent = experiment.get_with_parameters(db_session, uuid.uuid4())
    assert non_existent is None


def test_get_with_molecules(db_session, create_test_experiment, create_test_molecule):
    """Test retrieving an experiment with its associated molecules"""
    # Create a test experiment
    exp_type = ExperimentType(name="Molecule Test", description="Test molecules", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    test_exp = create_test_experiment(
        name="Test Get Molecules",
        created_by=1,  # Assuming user ID 1 exists
        experiment_type_id=exp_type.id
    )
    
    # Create test molecules and add them to the experiment
    mol1 = create_test_molecule(smiles="CC", created_by=1)
    mol2 = create_test_molecule(smiles="CCC", created_by=1)
    
    experiment.add_molecules(db_session, test_exp.id, [mol1.id, mol2.id])
    
    # Retrieve the experiment with molecules using get_with_molecules function
    retrieved_exp = experiment.get_with_molecules(db_session, test_exp.id)
    
    # Assert that the experiment has the correct number of molecules
    assert retrieved_exp.id == test_exp.id
    assert hasattr(retrieved_exp, "molecules")
    assert len(retrieved_exp.molecules) == 2
    
    # Assert that the molecules match the expected molecules
    molecule_smiles = {m.smiles for m in retrieved_exp.molecules}
    assert "CC" in molecule_smiles
    assert "CCC" in molecule_smiles
    
    # Test with a non-existent experiment ID
    non_existent = experiment.get_with_molecules(db_session, uuid.uuid4())
    assert non_existent is None


def test_get_by_status(db_session, create_test_experiment):
    """Test retrieving experiments with a specific status"""
    # Create multiple test experiments with different statuses
    exp_type = ExperimentType(name="Status Test", description="Test status", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    draft_exp1 = create_test_experiment(
        name="Draft Experiment 1",
        created_by=1,
        experiment_type_id=exp_type.id,
        status=ExperimentStatus.DRAFT.name
    )
    
    draft_exp2 = create_test_experiment(
        name="Draft Experiment 2",
        created_by=1,
        experiment_type_id=exp_type.id,
        status=ExperimentStatus.DRAFT.name
    )
    
    queued_exp = create_test_experiment(
        name="Queued Experiment",
        created_by=1,
        experiment_type_id=exp_type.id,
        status=ExperimentStatus.QUEUED.name
    )
    
    # Test retrieving experiments with a specific status
    draft_exps = experiment.get_by_status(db_session, ExperimentStatus.DRAFT.name)
    
    assert len(draft_exps) == 2
    draft_names = {exp.name for exp in draft_exps}
    assert "Draft Experiment 1" in draft_names
    assert "Draft Experiment 2" in draft_names
    
    # Test with pagination
    draft_exps_paged = experiment.get_by_status(db_session, ExperimentStatus.DRAFT.name, skip=1, limit=1)
    assert len(draft_exps_paged) == 1
    
    # Test with a status that no experiments have
    submitted_exps = experiment.get_by_status(db_session, ExperimentStatus.SUBMITTED.name)
    assert len(submitted_exps) == 0


def test_count_by_status(db_session, create_test_experiment):
    """Test counting experiments with a specific status"""
    # Create multiple test experiments with different statuses
    exp_type = ExperimentType(name="Count Status Test", description="Test count status", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    create_test_experiment(
        name="Draft Experiment Count 1",
        created_by=1,
        experiment_type_id=exp_type.id,
        status=ExperimentStatus.DRAFT.name
    )
    
    create_test_experiment(
        name="Draft Experiment Count 2",
        created_by=1,
        experiment_type_id=exp_type.id,
        status=ExperimentStatus.DRAFT.name
    )
    
    create_test_experiment(
        name="Queued Experiment Count",
        created_by=1,
        experiment_type_id=exp_type.id,
        status=ExperimentStatus.QUEUED.name
    )
    
    # Test counting experiments with a specific status
    draft_count = experiment.count_by_status(db_session, ExperimentStatus.DRAFT.name)
    assert draft_count == 2
    
    # Test with a status that no experiments have
    submitted_count = experiment.count_by_status(db_session, ExperimentStatus.SUBMITTED.name)
    assert submitted_count == 0


def test_update_status(db_session, create_test_experiment):
    """Test updating the status of an experiment"""
    # Create a test experiment with an initial status
    exp_type = ExperimentType(name="Update Status Test", description="Test update status", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    test_exp = create_test_experiment(
        name="Test Update Status",
        created_by=1,
        experiment_type_id=exp_type.id,
        status=ExperimentStatus.DRAFT.name
    )
    
    # Update the status using update_status function
    updated_exp = experiment.update_status(db_session, test_exp.id, ExperimentStatus.QUEUED.name)
    
    # Assert that the status was updated correctly
    assert updated_exp.id == test_exp.id
    assert updated_exp.status == ExperimentStatus.QUEUED
    
    # Test with a non-existent experiment ID
    non_existent = experiment.update_status(db_session, uuid.uuid4(), ExperimentStatus.QUEUED.name)
    assert non_existent is None


def test_get_by_user(db_session, create_test_experiment, test_pharma_user):
    """Test retrieving experiments created by a specific user"""
    # Create multiple test experiments with different user IDs
    exp_type = ExperimentType(name="User Test", description="Test user", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    user1_exp1 = create_test_experiment(
        name="User 1 Experiment 1",
        created_by=test_pharma_user.id,
        experiment_type_id=exp_type.id
    )
    
    user1_exp2 = create_test_experiment(
        name="User 1 Experiment 2",
        created_by=test_pharma_user.id,
        experiment_type_id=exp_type.id
    )
    
    user2_exp = create_test_experiment(
        name="User 2 Experiment",
        created_by=2,  # Assuming another user with ID 2
        experiment_type_id=exp_type.id
    )
    
    # Test retrieving experiments for a specific user
    user1_exps = experiment.get_by_user(db_session, test_pharma_user.id)
    
    assert len(user1_exps) == 2
    user1_names = {exp.name for exp in user1_exps}
    assert "User 1 Experiment 1" in user1_names
    assert "User 1 Experiment 2" in user1_names
    
    # Test with pagination
    user1_exps_paged = experiment.get_by_user(db_session, test_pharma_user.id, skip=1, limit=1)
    assert len(user1_exps_paged) == 1
    
    # Test with a user ID that has no experiments
    user3_exps = experiment.get_by_user(db_session, 3)  # Assuming no user with ID 3 or no experiments
    assert len(user3_exps) == 0


def test_count_by_user(db_session, create_test_experiment, test_pharma_user):
    """Test counting experiments created by a specific user"""
    # Create multiple test experiments with different user IDs
    exp_type = ExperimentType(name="Count User Test", description="Test count user", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    create_test_experiment(
        name="User 1 Count Experiment 1",
        created_by=test_pharma_user.id,
        experiment_type_id=exp_type.id
    )
    
    create_test_experiment(
        name="User 1 Count Experiment 2",
        created_by=test_pharma_user.id,
        experiment_type_id=exp_type.id
    )
    
    create_test_experiment(
        name="User 2 Count Experiment",
        created_by=2,  # Assuming another user with ID 2
        experiment_type_id=exp_type.id
    )
    
    # Test counting experiments for a specific user
    user1_count = experiment.count_by_user(db_session, test_pharma_user.id)
    assert user1_count == 2
    
    # Test with a user ID that has no experiments
    user3_count = experiment.count_by_user(db_session, 3)  # Assuming no user with ID 3 or no experiments
    assert user3_count == 0


def test_add_molecules(db_session, create_test_experiment, create_test_molecule):
    """Test adding molecules to an experiment"""
    # Create a test experiment
    exp_type = ExperimentType(name="Add Molecules Test", description="Test add molecules", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    test_exp = create_test_experiment(
        name="Test Add Molecules",
        created_by=1,
        experiment_type_id=exp_type.id
    )
    
    # Create test molecules
    mol1 = create_test_molecule(smiles="CCO", created_by=1)
    mol2 = create_test_molecule(smiles="CCCO", created_by=1)
    
    # Add molecules to the experiment using add_molecules function
    result = experiment.add_molecules(db_session, test_exp.id, [mol1.id, mol2.id])
    
    # Assert that the molecules were added successfully
    assert result["success"] is True
    assert result["added"] == 2
    assert result["failed"] == 0
    
    # Test adding molecules that are already in the experiment
    result2 = experiment.add_molecules(db_session, test_exp.id, [mol1.id])
    assert result2["success"] is True
    assert result2["added"] == 1  # Counted as success but not actually added again
    
    # Test with a non-existent experiment ID
    result3 = experiment.add_molecules(db_session, uuid.uuid4(), [mol1.id])
    assert result3["success"] is False
    assert "Experiment not found" in result3["error"]


def test_remove_molecules(db_session, create_test_experiment, create_test_molecule):
    """Test removing molecules from an experiment"""
    # Create a test experiment
    exp_type = ExperimentType(name="Remove Molecules Test", description="Test remove molecules", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    test_exp = create_test_experiment(
        name="Test Remove Molecules",
        created_by=1,
        experiment_type_id=exp_type.id
    )
    
    # Create test molecules and add them to the experiment
    mol1 = create_test_molecule(smiles="CCO", created_by=1)
    mol2 = create_test_molecule(smiles="CCCO", created_by=1)
    mol3 = create_test_molecule(smiles="CCCCO", created_by=1)
    
    experiment.add_molecules(db_session, test_exp.id, [mol1.id, mol2.id, mol3.id])
    
    # Remove molecules from the experiment using remove_molecules function
    result = experiment.remove_molecules(db_session, test_exp.id, [mol1.id, mol2.id])
    
    # Assert that the molecules were removed successfully
    assert result["success"] is True
    assert result["removed"] == 2
    assert result["failed"] == 0
    
    # Test removing molecules that are not in the experiment
    result2 = experiment.remove_molecules(db_session, test_exp.id, [mol1.id])
    assert result2["success"] is True
    assert result2["removed"] == 0
    assert result2["failed"] == 1
    
    # Test with a non-existent experiment ID
    result3 = experiment.remove_molecules(db_session, uuid.uuid4(), [mol3.id])
    assert result3["success"] is False
    assert "Experiment not found" in result3["error"]


def test_get_molecules(db_session, create_test_experiment, create_test_molecule):
    """Test retrieving molecules associated with an experiment"""
    # Create a test experiment
    exp_type = ExperimentType(name="Get Molecules Test", description="Test get molecules", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    test_exp = create_test_experiment(
        name="Test Get Molecules Func",
        created_by=1,
        experiment_type_id=exp_type.id
    )
    
    # Create test molecules and add them to the experiment
    molecules = []
    for i in range(5):
        mol = create_test_molecule(smiles=f"C{'C' * i}O", created_by=1)
        molecules.append(mol)
    
    experiment.add_molecules(db_session, test_exp.id, [m.id for m in molecules])
    
    # Retrieve molecules using get_molecules function
    exp_molecules = experiment.get_molecules(db_session, test_exp.id)
    assert len(exp_molecules) == 5
    
    # Test with pagination
    paged_molecules = experiment.get_molecules(db_session, test_exp.id, skip=2, limit=2)
    assert len(paged_molecules) == 2
    
    # Test with a non-existent experiment ID
    non_existent = experiment.get_molecules(db_session, uuid.uuid4())
    assert len(non_existent) == 0


def test_count_molecules(db_session, create_test_experiment, create_test_molecule):
    """Test counting molecules associated with an experiment"""
    # Create a test experiment
    exp_type = ExperimentType(name="Count Molecules Test", description="Test count molecules", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    test_exp = create_test_experiment(
        name="Test Count Molecules",
        created_by=1,
        experiment_type_id=exp_type.id
    )
    
    # Create test molecules and add them to the experiment
    molecules = []
    for i in range(3):
        mol = create_test_molecule(smiles=f"C{'C' * i}N", created_by=1)
        molecules.append(mol)
    
    experiment.add_molecules(db_session, test_exp.id, [m.id for m in molecules])
    
    # Count molecules using count_molecules function
    count = experiment.count_molecules(db_session, test_exp.id)
    assert count == 3
    
    # Test with a non-existent experiment ID
    non_existent_count = experiment.count_molecules(db_session, uuid.uuid4())
    assert non_existent_count == 0


def test_is_molecule_in_experiment(db_session, create_test_experiment, create_test_molecule):
    """Test checking if a molecule is associated with an experiment"""
    # Create a test experiment
    exp_type = ExperimentType(name="Molecule In Exp Test", description="Test molecule in exp", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    test_exp = create_test_experiment(
        name="Test Is Molecule In Experiment",
        created_by=1,
        experiment_type_id=exp_type.id
    )
    
    # Create test molecules and add some to the experiment
    mol1 = create_test_molecule(smiles="CCBr", created_by=1)
    mol2 = create_test_molecule(smiles="CCCl", created_by=1)
    mol3 = create_test_molecule(smiles="CCF", created_by=1)
    
    experiment.add_molecules(db_session, test_exp.id, [mol1.id, mol2.id])
    
    # Check if molecules are in the experiment using is_molecule_in_experiment function
    assert experiment.is_molecule_in_experiment(db_session, test_exp.id, mol1.id) is True
    assert experiment.is_molecule_in_experiment(db_session, test_exp.id, mol2.id) is True
    assert experiment.is_molecule_in_experiment(db_session, test_exp.id, mol3.id) is False
    
    # Test with non-existent experiment or molecule IDs
    assert experiment.is_molecule_in_experiment(db_session, uuid.uuid4(), mol1.id) is False
    assert experiment.is_molecule_in_experiment(db_session, test_exp.id, uuid.uuid4()) is False


def test_get_by_type(db_session, create_test_experiment):
    """Test retrieving experiments of a specific type"""
    # Create multiple test experiments with different types
    type1 = ExperimentType(name="Type 1", description="Test type 1", category="Test")
    type2 = ExperimentType(name="Type 2", description="Test type 2", category="Test")
    db_session.add_all([type1, type2])
    db_session.flush()
    
    type1_exp1 = create_test_experiment(
        name="Type 1 Experiment 1",
        created_by=1,
        experiment_type_id=type1.id
    )
    
    type1_exp2 = create_test_experiment(
        name="Type 1 Experiment 2",
        created_by=1,
        experiment_type_id=type1.id
    )
    
    type2_exp = create_test_experiment(
        name="Type 2 Experiment",
        created_by=1,
        experiment_type_id=type2.id
    )
    
    # Test retrieving experiments with a specific type
    type1_exps = experiment.get_by_type(db_session, type1.id)
    
    assert len(type1_exps) == 2
    type1_names = {exp.name for exp in type1_exps}
    assert "Type 1 Experiment 1" in type1_names
    assert "Type 1 Experiment 2" in type1_names
    
    # Test with pagination
    type1_exps_paged = experiment.get_by_type(db_session, type1.id, skip=1, limit=1)
    assert len(type1_exps_paged) == 1
    
    # Test with a type that no experiments have
    type3_exps = experiment.get_by_type(db_session, uuid.uuid4())
    assert len(type3_exps) == 0


def test_get_parameter_value(db_session, create_test_experiment):
    """Test retrieving the value of a specific parameter for an experiment"""
    # Create a test experiment with parameters
    exp_type = ExperimentType(name="Param Value Test", description="Test param value", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    test_exp = create_test_experiment(
        name="Test Get Parameter Value",
        created_by=1,
        experiment_type_id=exp_type.id,
        parameters={
            "target": "Protein X",
            "concentration": "5 µM"
        }
    )
    
    # Retrieve parameter value using get_parameter_value function
    target_value = experiment.get_parameter_value(db_session, test_exp.id, "target")
    concentration_value = experiment.get_parameter_value(db_session, test_exp.id, "concentration")
    
    # Assert that the correct parameter value is returned
    assert target_value == "Protein X"
    assert concentration_value == "5 µM"
    
    # Test with a non-existent parameter name
    non_existent_param = experiment.get_parameter_value(db_session, test_exp.id, "non_existent")
    assert non_existent_param is None
    
    # Test with a non-existent experiment ID
    non_existent_exp = experiment.get_parameter_value(db_session, uuid.uuid4(), "target")
    assert non_existent_exp is None


def test_set_parameter_value(db_session, create_test_experiment):
    """Test setting the value of a specific parameter for an experiment"""
    # Create a test experiment with parameters
    exp_type = ExperimentType(name="Set Param Test", description="Test set param", category="Test")
    db_session.add(exp_type)
    db_session.flush()
    
    test_exp = create_test_experiment(
        name="Test Set Parameter Value",
        created_by=1,
        experiment_type_id=exp_type.id,
        parameters={
            "temperature": "25°C"
        }
    )
    
    # Set parameter value using set_parameter_value function
    updated_param = experiment.set_parameter_value(db_session, test_exp.id, "temperature", "30°C")
    
    # Assert that the parameter value was updated correctly
    assert updated_param.parameter_name == "temperature"
    assert updated_param.parameter_value == "30°C"
    
    # Test setting a new parameter that doesn't exist yet
    new_param = experiment.set_parameter_value(db_session, test_exp.id, "pH", "7.4")
    
    # Assert that the new parameter was created correctly
    assert new_param.parameter_name == "pH"
    assert new_param.parameter_value == "7.4"
    
    # Test with a non-existent experiment ID
    with pytest.raises(Exception):  # Should raise some exception
        experiment.set_parameter_value(db_session, uuid.uuid4(), "temperature", "35°C")


def test_filter_experiments(db_session, create_test_experiment, create_test_molecule):
    """Test filtering experiments based on various criteria"""
    # Create multiple test experiments with different attributes
    binding_type = ExperimentType(name="Binding", description="Binding assay", category="Assay")
    adme_type = ExperimentType(name="ADME", description="ADME assessment", category="Assay")
    db_session.add_all([binding_type, adme_type])
    db_session.flush()
    
    # Create molecules
    mol1 = create_test_molecule(smiles="CCO", created_by=1)
    mol2 = create_test_molecule(smiles="c1ccccc1", created_by=1)
    
    # Create experiments with various attributes
    exp1 = create_test_experiment(
        name="Binding Study 1",
        created_by=1,
        experiment_type_id=binding_type.id,
        status=ExperimentStatus.DRAFT.name
    )
    
    exp2 = create_test_experiment(
        name="Binding Study 2",
        created_by=1,
        experiment_type_id=binding_type.id,
        status=ExperimentStatus.QUEUED.name
    )
    
    exp3 = create_test_experiment(
        name="ADME Assessment",
        created_by=2,
        experiment_type_id=adme_type.id,
        status=ExperimentStatus.DRAFT.name
    )
    
    # Add molecules to experiments
    experiment.add_molecules(db_session, exp1.id, [mol1.id])
    experiment.add_molecules(db_session, exp2.id, [mol1.id, mol2.id])
    experiment.add_molecules(db_session, exp3.id, [mol2.id])
    
    # Test filtering by name
    name_filter = {"name": "Binding"}
    binding_exps, binding_count = experiment.filter_experiments(db_session, name_filter)
    assert binding_count == 2
    binding_names = {exp.name for exp in binding_exps}
    assert "Binding Study 1" in binding_names
    assert "Binding Study 2" in binding_names
    
    # Test filtering by status
    status_filter = {"status": ExperimentStatus.DRAFT.name}
    draft_exps, draft_count = experiment.filter_experiments(db_session, status_filter)
    assert draft_count == 2
    draft_ids = {exp.id for exp in draft_exps}
    assert exp1.id in draft_ids
    assert exp3.id in draft_ids
    
    # Test filtering by type
    type_filter = {"type_id": binding_type.id}
    binding_type_exps, binding_type_count = experiment.filter_experiments(db_session, type_filter)
    assert binding_type_count == 2
    binding_type_ids = {exp.id for exp in binding_type_exps}
    assert exp1.id in binding_type_ids
    assert exp2.id in binding_type_ids
    
    # Test filtering by user
    user_filter = {"created_by": 2}
    user_exps, user_count = experiment.filter_experiments(db_session, user_filter)
    assert user_count == 1
    assert user_exps[0].id == exp3.id
    
    # Test filtering by molecule
    molecule_filter = {"molecule_id": mol2.id}
    mol_exps, mol_count = experiment.filter_experiments(db_session, molecule_filter)
    assert mol_count == 2
    mol_ids = {exp.id for exp in mol_exps}
    assert exp2.id in mol_ids
    assert exp3.id in mol_ids
    
    # Test filtering with multiple criteria
    complex_filter = {"status": ExperimentStatus.DRAFT.name, "type_id": binding_type.id}
    complex_exps, complex_count = experiment.filter_experiments(db_session, complex_filter)
    assert complex_count == 1
    assert complex_exps[0].id == exp1.id
    
    # Test with sorting (ascending and descending)
    sort_filter_desc = {"sort_by": "name", "sort_desc": True}
    sorted_desc_exps, _ = experiment.filter_experiments(db_session, sort_filter_desc)
    assert sorted_desc_exps[0].name > sorted_desc_exps[-1].name
    
    sort_filter_asc = {"sort_by": "name", "sort_desc": False}
    sorted_asc_exps, _ = experiment.filter_experiments(db_session, sort_filter_asc)
    assert sorted_asc_exps[0].name < sorted_asc_exps[-1].name
    
    # Test with pagination
    paged_filter = {}
    paged_exps, total = experiment.filter_experiments(db_session, paged_filter, skip=1, limit=1)
    assert len(paged_exps) == 1
    assert total == 3  # Total count should still be the total number of experiments