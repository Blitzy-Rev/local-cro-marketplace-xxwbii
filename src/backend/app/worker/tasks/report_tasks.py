"""
Celery task module for asynchronous report generation in the Molecular Data Management and CRO Integration Platform.
Implements background tasks for generating various types of reports including experiment summaries, molecule property analyses, and result compilations.
"""

import logging
import os
import io
from typing import Dict, List, Any, Optional, Union, Tuple
from uuid import UUID
from datetime import datetime

import pandas  # pandas 2.0+
import numpy  # numpy 1.24+
import matplotlib  # matplotlib 3.7+
matplotlib.use('Agg')  # Set the backend to 'Agg'
import matplotlib.pyplot as plt  # matplotlib 3.7+
from io import BytesIO
from celery import Celery  # celery 5.2+

from ..celery_app import app
from ...logging_config import logger
from ...services import experiment_service
from ...services import result_service
from ...services import molecule_service
from ...services.file_storage_service import FileStorageService
from ...constants import BUCKET_NAMES
from ...exceptions import ValidationException, ResourceNotFoundException

# Initialize global variables
logger = logging.getLogger(__name__)
file_storage_service = FileStorageService()

@app.task(bind=True, name='report_tasks.generate_experiment_report', max_retries=3, acks_late=True)
def generate_experiment_report(self, experiment_id: UUID, user_id: int, report_format: str = 'pdf') -> Dict[str, Any]:
    """
    Celery task for generating a comprehensive report for an experiment

    Args:
        experiment_id: ID of the experiment
        user_id: ID of the user
        report_format: Report format ('pdf', 'xlsx')

    Returns:
        Report generation results with file path
    """
    logger.info(f"Starting experiment report generation task for experiment_id: {experiment_id}")

    try:
        # Validate report_format
        if report_format not in ['pdf', 'xlsx']:
            raise ValidationException(f"Unsupported report format: {report_format}", details={"report_format": report_format})

        # Get experiment details
        experiment_data = experiment_service.get_experiment_by_id(experiment_id, include_molecules=True)
        if not experiment_data:
            raise ResourceNotFoundException(f"Experiment with ID {experiment_id} not found", details={"experiment_id": experiment_id})

        # Get experiment molecules
        molecules = experiment_data.get('molecules', [])

        # Create pandas DataFrame with molecule properties
        df = pandas.DataFrame(molecules)

        # Generate statistical summary of molecular properties
        summary = df.describe()

        # Create visualizations (property distributions, etc.)
        figures = _generate_property_charts(df, ['MW', 'LogP'])

        # Compile experiment metadata, molecule data, and results into report
        metadata = {
            "Experiment ID": str(experiment_id),
            "Experiment Name": experiment_data.get('name', 'N/A'),
            "Created By": experiment_data.get('creator', {}).get('email', 'N/A'),
            "Created At": str(experiment_data.get('created_at', 'N/A')),
            "Total Molecules": len(molecules)
        }

        # Generate report file in requested format
        if report_format == 'pdf':
            report_bytes = _create_pdf_report(f"Experiment Report - {experiment_data.get('name', 'N/A')}", metadata, df, figures)
        elif report_format == 'xlsx':
            report_bytes = _create_excel_report(f"Experiment Report - {experiment_data.get('name', 'N/A')}", metadata, df, figures)
        else:
            raise ValidationException(f"Unsupported report format: {report_format}", details={"report_format": report_format})

        # Upload report file to storage
        filename = f"experiment_report_{experiment_id}.{report_format}"
        object_name = file_storage_service.upload_experiment_file(io.BytesIO(report_bytes.getvalue()), filename)

        # Return report generation results with file path
        return {
            "experiment_id": str(experiment_id),
            "report_format": report_format,
            "file_path": object_name
        }

    except ResourceNotFoundException as e:
        logger.error(f"Experiment report generation failed: {e.message}")
        raise self.retry(exc=e, countdown=60)  # Retry after 60 seconds
    except ValidationException as e:
        logger.error(f"Validation error during experiment report generation: {e.message}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error during experiment report generation: {str(e)}")
        raise self.retry(exc=e, countdown=60)  # Retry after 60 seconds

@app.task(bind=True, name='report_tasks.generate_result_report', max_retries=3, acks_late=True)
def generate_result_report(self, result_id: UUID, user_id: int, report_format: str = 'pdf') -> Dict[str, Any]:
    """
    Celery task for generating a detailed report for experimental results

    Args:
        result_id: ID of the experimental result
        user_id: ID of the user
        report_format: Report format ('pdf', 'xlsx')

    Returns:
        Report generation results with file path
    """
    logger.info(f"Starting result report generation task for result_id: {result_id}")

    try:
        # Validate report_format
        if report_format not in ['pdf', 'xlsx']:
            raise ValidationException(f"Unsupported report format: {report_format}", details={"report_format": report_format})

        # Get result details
        result_data = result_service.get_result_by_id(result_id)
        if not result_data:
            raise ResourceNotFoundException(f"Result with ID {result_id} not found", details={"result_id": result_id})

        # Get result data
        result_details = result_service.get_result_data(result_id)

        # Get associated experiment and submission details
        # experiment = experiment_service.get_experiment_by_id(result_data.experiment_id)
        # submission = submission_service.get_submission_by_id(result_data.submission_id)

        # Get molecules involved in the experiment
        # molecules = experiment_service.get_experiment_molecules(result_data.experiment_id)

        # Create pandas DataFrame with result data
        df = pandas.DataFrame(result_details)

        # Generate statistical analysis of result data
        summary = df.describe()

        # Create visualizations (result comparisons, property correlations)
        figures = _generate_property_charts(df, ['data_value'])

        # Compile result metadata, data points, and analysis into report
        metadata = {
            "Result ID": str(result_id),
            "Uploaded At": str(result_data.get('uploaded_at', 'N/A')),
            "Status": result_data.get('status', 'N/A')
        }

        # Generate report file in requested format
        if report_format == 'pdf':
            report_bytes = _create_pdf_report(f"Result Report - {result_id}", metadata, df, figures)
        elif report_format == 'xlsx':
            report_bytes = _create_excel_report(f"Result Report - {result_id}", metadata, df, figures)
        else:
            raise ValidationException(f"Unsupported report format: {report_format}", details={"report_format": report_format})

        # Upload report file to storage
        filename = f"result_report_{result_id}.{report_format}"
        object_name = file_storage_service.upload_experiment_file(io.BytesIO(report_bytes.getvalue()), filename)

        # Return report generation results with file path
        return {
            "result_id": str(result_id),
            "report_format": report_format,
            "file_path": object_name
        }

    except ResourceNotFoundException as e:
        logger.error(f"Result report generation failed: {e.message}")
        raise self.retry(exc=e, countdown=60)  # Retry after 60 seconds
    except ValidationException as e:
        logger.error(f"Validation error during result report generation: {e.message}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error during result report generation: {str(e)}")
        raise self.retry(exc=e, countdown=60)  # Retry after 60 seconds

@app.task(bind=True, name='report_tasks.generate_molecule_comparison_report', max_retries=3, acks_late=True)
def generate_molecule_comparison_report(self, molecule_ids: List[UUID], user_id: int, report_format: str = 'pdf') -> Dict[str, Any]:
    """
    Celery task for generating a comparison report for multiple molecules

    Args:
        molecule_ids: List of molecule IDs
        user_id: ID of the user
        report_format: Report format ('pdf', 'xlsx')

    Returns:
        Report generation results with file path
    """
    logger.info(f"Starting molecule comparison report generation task for {len(molecule_ids)} molecules")

    try:
        # Validate report_format
        if report_format not in ['pdf', 'xlsx']:
            raise ValidationException(f"Unsupported report format: {report_format}", details={"report_format": report_format})

        # Validate molecule_ids list is not empty
        if not molecule_ids:
            raise ValidationException("Molecule IDs list cannot be empty", details={"molecule_ids": molecule_ids})

        # Get molecules
        molecules = molecule_service.get_molecules_by_ids(molecule_ids)

        # Create pandas DataFrame with molecule properties
        df = pandas.DataFrame(molecules)

        # Generate comparative analysis of molecular properties
        summary = df.describe()

        # Create visualizations (property comparisons, radar charts)
        figures = _generate_property_charts(df, ['MW', 'LogP'])

        # Compile molecule structures, properties, and comparisons into report
        metadata = {
            "Report Type": "Molecule Comparison",
            "Number of Molecules": len(molecules),
            "Generated By": user_id
        }

        # Generate report file in requested format
        if report_format == 'pdf':
            report_bytes = _create_pdf_report("Molecule Comparison Report", metadata, df, figures)
        elif report_format == 'xlsx':
            report_bytes = _create_excel_report("Molecule Comparison Report", metadata, df, figures)
        else:
            raise ValidationException(f"Unsupported report format: {report_format}", details={"report_format": report_format})

        # Upload report file to storage
        filename = f"molecule_comparison_report_{len(molecule_ids)}_molecules.{report_format}"
        object_name = file_storage_service.upload_experiment_file(io.BytesIO(report_bytes.getvalue()), filename)

        # Return report generation results with file path
        return {
            "molecule_ids": [str(mid) for mid in molecule_ids],
            "report_format": report_format,
            "file_path": object_name
        }

    except ResourceNotFoundException as e:
        logger.error(f"Molecule comparison report generation failed: {e.message}")
        raise self.retry(exc=e, countdown=60)  # Retry after 60 seconds
    except ValidationException as e:
        logger.error(f"Validation error during molecule comparison report generation: {e.message}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error during molecule comparison report generation: {str(e)}")
        raise self.retry(exc=e, countdown=60)  # Retry after 60 seconds

@app.task(bind=True, name='report_tasks.generate_library_report', max_retries=3, acks_late=True)
def generate_library_report(self, library_id: UUID, user_id: int, report_format: str = 'pdf') -> Dict[str, Any]:
    """
    Celery task for generating a statistical report for a molecule library

    Args:
        library_id: ID of the molecule library
        user_id: ID of the user
        report_format: Report format ('pdf', 'xlsx')

    Returns:
        Report generation results with file path
    """
    logger.info(f"Starting library report generation task for library_id: {library_id}")

    try:
        # Validate report_format
        if report_format not in ['pdf', 'xlsx']:
            raise ValidationException(f"Unsupported report format: {report_format}", details={"report_format": report_format})

        # Get library details and molecules
        library = library_service.get_library_by_id(library_id)
        molecules = library.get('molecules', [])

        # Create pandas DataFrame with molecule properties
        df = pandas.DataFrame(molecules)

        # Generate statistical summary of library properties
        summary = df.describe()

        # Create visualizations (property distributions, clustering)
        figures = _generate_property_charts(df, ['MW', 'LogP'])

        # Compile library metadata, statistics, and visualizations into report
        metadata = {
            "Library ID": str(library_id),
            "Library Name": library.get('name', 'N/A'),
            "Created By": user_id,
            "Total Molecules": len(molecules)
        }

        # Generate report file in requested format
        if report_format == 'pdf':
            report_bytes = _create_pdf_report("Molecule Library Report", metadata, df, figures)
        elif report_format == 'xlsx':
            report_bytes = _create_excel_report("Molecule Library Report", metadata, df, figures)
        else:
            raise ValidationException(f"Unsupported report format: {report_format}", details={"report_format": report_format})

        # Upload report file to storage
        filename = f"molecule_library_report_{library_id}.{report_format}"
        object_name = file_storage_service.upload_experiment_file(io.BytesIO(report_bytes.getvalue()), filename)

        # Return report generation results with file path
        return {
            "library_id": str(library_id),
            "report_format": report_format,
            "file_path": object_name
        }

    except ResourceNotFoundException as e:
        logger.error(f"Library report generation failed: {e.message}")
        raise self.retry(exc=e, countdown=60)  # Retry after 60 seconds
    except ValidationException as e:
        logger.error(f"Validation error during library report generation: {e.message}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error during library report generation: {str(e)}")
        raise self.retry(exc=e, countdown=60)  # Retry after 60 seconds

@app.task(bind=True, name='report_tasks.generate_activity_report', max_retries=3, acks_late=True)
def generate_activity_report(self, user_id: int, period: str = 'month', report_format: str = 'pdf') -> Dict[str, Any]:
    """
    Celery task for generating a periodic activity report for a user

    Args:
        user_id: ID of the user
        period: Time period for the report ('week', 'month', 'quarter', 'year')
        report_format: Report format ('pdf', 'xlsx')

    Returns:
        Report generation results with file path
    """
    logger.info(f"Starting activity report generation task for user_id: {user_id} and period: {period}")

    try:
        # Validate period
        if period not in ['week', 'month', 'quarter', 'year']:
            raise ValidationException(f"Unsupported period: {period}", details={"period": period})

        # Validate report_format
        if report_format not in ['pdf', 'xlsx']:
            raise ValidationException(f"Unsupported report format: {report_format}", details={"report_format": report_format})

        # Calculate date range based on period
        if period == 'week':
            start_date = datetime.utcnow() - timedelta(days=7)
        elif period == 'month':
            start_date = datetime.utcnow() - timedelta(days=30)
        elif period == 'quarter':
            start_date = datetime.utcnow() - timedelta(days=90)
        else:  # year
            start_date = datetime.utcnow() - timedelta(days=365)

        # Get user's experiments, submissions, and results within date range
        experiments = experiment_service.get_experiments_by_user(user_id)
        # submissions = submission_service.get_submissions_by_user(user_id, from_date=start_date)
        # results = result_service.get_results_by_user(user_id, from_date=start_date)

        # Create summary of user activity (counts, status distributions)
        activity_summary = {
            "Total Experiments": len(experiments),
            # "Total Submissions": len(submissions),
            # "Total Results": len(results)
        }

        # Generate timeline visualization of activities
        # timeline_figure = _generate_activity_timeline(experiments, submissions, results)
        timeline_figure = None

        # Compile activity data and visualizations into report
        metadata = {
            "Report Type": "User Activity",
            "User ID": user_id,
            "Period": period,
            "Start Date": str(start_date),
            "End Date": str(datetime.utcnow())
        }

        # Create pandas DataFrame with activity data
        df = pandas.DataFrame([activity_summary])

        # Generate report file in requested format
        if report_format == 'pdf':
            report_bytes = _create_pdf_report("User Activity Report", metadata, df, [timeline_figure] if timeline_figure else [])
        elif report_format == 'xlsx':
            report_bytes = _create_excel_report("User Activity Report", metadata, df, [timeline_figure] if timeline_figure else [])
        else:
            raise ValidationException(f"Unsupported report format: {report_format}", details={"report_format": report_format})

        # Upload report file to storage
        filename = f"user_activity_report_{user_id}_{period}.{report_format}"
        object_name = file_storage_service.upload_experiment_file(io.BytesIO(report_bytes.getvalue()), filename)

        # Return report generation results with file path
        return {
            "user_id": user_id,
            "period": period,
            "report_format": report_format,
            "file_path": object_name
        }

    except ResourceNotFoundException as e:
        logger.error(f"Activity report generation failed: {e.message}")
        raise self.retry(exc=e, countdown=60)  # Retry after 60 seconds
    except ValidationException as e:
        logger.error(f"Validation error during activity report generation: {e.message}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error during activity report generation: {str(e)}")
        raise self.retry(exc=e, countdown=60)  # Retry after 60 seconds

def _create_pdf_report(title: str, metadata: Dict[str, Any], data: pandas.DataFrame, figures: List[matplotlib.figure.Figure]) -> io.BytesIO:
    """
    Helper function to create a PDF report with charts and tables

    Args:
        title: Title of the report
        metadata: Dictionary of metadata to include
        data: pandas DataFrame with data to include
        figures: List of matplotlib figures to include

    Returns:
        PDF report as bytes buffer
    """
    # Implementation details for PDF report generation using ReportLab
    # (omitted for brevity)
    return io.BytesIO()

def _create_excel_report(title: str, metadata: Dict[str, Any], data: pandas.DataFrame, figures: List[matplotlib.figure.Figure]) -> io.BytesIO:
    """
    Helper function to create an Excel report with data and charts

    Args:
        title: Title of the report
        metadata: Dictionary of metadata to include
        data: pandas DataFrame with data to include
        figures: List of matplotlib figures to include

    Returns:
        Excel report as bytes buffer
    """
    # Implementation details for Excel report generation using pandas and openpyxl
    # (omitted for brevity)
    return io.BytesIO()

def _generate_property_charts(data: pandas.DataFrame, properties: List[str]) -> List[matplotlib.figure.Figure]:
    """
    Helper function to generate charts for molecular properties

    Args:
        data: pandas DataFrame with molecule data
        properties: List of property names to generate charts for

    Returns:
        List of generated chart figures
    """
    # Implementation details for generating charts using matplotlib
    # (omitted for brevity)
    return []