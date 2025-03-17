from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy import text
import redis
import minio
import sqlalchemy.exc
import time
from typing import Dict, Any

from ..deps import get_db_session
from ..../../db.session import get_engine
from ..../../core.config import get_redis_settings, get_minio_settings
from ..../../worker.celery_app import app
from ..../../logging_config import logger

# Create a router for health check endpoints
health_router = APIRouter(prefix='/health', tags=['health'])


def check_database() -> Dict[str, Any]:
    """
    Checks database connectivity by executing a simple query
    
    Returns:
        Dict[str, Any]: Database health status with connection state and latency
    """
    try:
        # Get database engine
        engine = get_engine()
        
        # Measure query execution time
        start_time = time.time()
        # Execute a simple query to verify connection
        with engine.connect() as conn:
            conn.execute(text("SELECT 1"))
            conn.commit()
        end_time = time.time()
        
        # Calculate query latency in milliseconds
        latency_ms = round((end_time - start_time) * 1000, 2)
        
        # Return success status with latency information
        return {
            "status": "UP",
            "latency_ms": latency_ms
        }
    except sqlalchemy.exc.SQLAlchemyError as e:
        # Log database connection error
        logger.error(f"Database health check failed: {str(e)}")
        
        # Return error status with details
        return {
            "status": "DOWN",
            "error": str(e)
        }


def check_redis() -> Dict[str, Any]:
    """
    Checks Redis connectivity by performing a ping operation
    
    Returns:
        Dict[str, Any]: Redis health status with connection state and latency
    """
    try:
        # Get Redis settings
        redis_settings = get_redis_settings()
        
        # Create Redis client
        redis_client = redis.Redis(
            host=redis_settings["host"],
            port=redis_settings["port"],
            password=redis_settings["password"],
            db=redis_settings["db"],
            socket_timeout=5
        )
        
        # Measure ping operation time
        start_time = time.time()
        redis_client.ping()
        end_time = time.time()
        
        # Calculate operation latency in milliseconds
        latency_ms = round((end_time - start_time) * 1000, 2)
        
        # Return success status with latency information
        return {
            "status": "UP",
            "latency_ms": latency_ms
        }
    except redis.RedisError as e:
        # Log Redis connection error
        logger.error(f"Redis health check failed: {str(e)}")
        
        # Return error status with details
        return {
            "status": "DOWN",
            "error": str(e)
        }


def check_minio() -> Dict[str, Any]:
    """
    Checks MinIO connectivity by listing buckets
    
    Returns:
        Dict[str, Any]: MinIO health status with connection state and latency
    """
    try:
        # Get MinIO settings
        minio_settings = get_minio_settings()
        
        # Create MinIO client
        minio_client = minio.Minio(
            f"{minio_settings['host']}:{minio_settings['port']}",
            access_key=minio_settings["access_key"],
            secret_key=minio_settings["secret_key"],
            secure=minio_settings["secure"]
        )
        
        # Measure list_buckets operation time
        start_time = time.time()
        minio_client.list_buckets()
        end_time = time.time()
        
        # Calculate operation latency in milliseconds
        latency_ms = round((end_time - start_time) * 1000, 2)
        
        # Return success status with latency information
        return {
            "status": "UP",
            "latency_ms": latency_ms
        }
    except Exception as e:
        # Log MinIO connection error
        logger.error(f"MinIO health check failed: {str(e)}")
        
        # Return error status with details
        return {
            "status": "DOWN",
            "error": str(e)
        }


def check_celery() -> Dict[str, Any]:
    """
    Checks Celery worker availability by inspecting active workers
    
    Returns:
        Dict[str, Any]: Celery health status with worker availability
    """
    try:
        # Create inspector from Celery app
        inspector = app.control.inspect()
        
        # Get active workers
        active_workers = inspector.active()
        
        # Check if any workers are active
        if active_workers and len(active_workers) > 0:
            # Return success status with worker count
            return {
                "status": "UP",
                "workers": len(active_workers)
            }
        else:
            # No active workers
            return {
                "status": "DOWN",
                "error": "No active workers found"
            }
    except Exception as e:
        # Log Celery connection error
        logger.error(f"Celery health check failed: {str(e)}")
        
        # Return error status with details
        return {
            "status": "DOWN",
            "error": str(e)
        }


def get_system_metrics() -> Dict[str, Any]:
    """
    Collects basic system metrics like memory and CPU usage
    
    Returns:
        Dict[str, Any]: System metrics including memory and CPU usage
    """
    try:
        # Try to import psutil if available
        import psutil
        
        # Collect memory usage information
        memory = psutil.virtual_memory()
        memory_metrics = {
            "total": memory.total,
            "used": memory.used,
            "percent": memory.percent
        }
        
        # Collect CPU usage information
        cpu_metrics = {
            "percent": psutil.cpu_percent(interval=0.1)
        }
        
        # Collect disk usage information
        disk = psutil.disk_usage('/')
        disk_metrics = {
            "total": disk.total,
            "used": disk.used,
            "percent": disk.percent
        }
        
        # Return all metrics
        return {
            "memory": memory_metrics,
            "cpu": cpu_metrics,
            "disk": disk_metrics
        }
    except ImportError:
        # psutil not available, return limited metrics
        logger.warning("psutil not installed, returning limited system metrics")
        return {
            "memory": {"total": 0, "used": 0, "percent": 0},
            "cpu": {"percent": 0},
            "disk": {"total": 0, "used": 0, "percent": 0}
        }


@health_router.get("/live", summary="Simple liveness check")
async def liveness_check():
    """
    Simple liveness check that returns UP status if the service is running.
    This is a lightweight check that doesn't verify any dependencies.
    
    Returns:
        dict: Simple status response
    """
    return {"status": "UP"}


@health_router.get("/ready", summary="Readiness check")
async def readiness_check():
    """
    Readiness check that verifies if the service is ready to handle requests.
    Checks database and Redis connectivity as these are essential for operation.
    
    Returns:
        dict: Readiness status with details about each dependency
    """
    # Check database and Redis
    db_status = check_database()
    redis_status = check_redis()
    
    # Determine overall status (UP only if all dependencies are UP)
    overall_status = "UP" if db_status["status"] == "UP" and redis_status["status"] == "UP" else "DOWN"
    
    # Return combined status
    return {
        "status": overall_status,
        "details": {
            "database": db_status,
            "redis": redis_status
        }
    }


@health_router.get("/deep", summary="Deep health check")
async def deep_health_check():
    """
    Comprehensive health check of all system components including database,
    Redis, MinIO, and Celery workers.
    
    Returns:
        dict: Detailed health status of all components
    """
    # Check all components
    db_status = check_database()
    redis_status = check_redis()
    minio_status = check_minio()
    celery_status = check_celery()
    
    # Determine overall status (UP only if all components are UP)
    overall_status = "UP"
    if any(component["status"] == "DOWN" for component in [db_status, redis_status, minio_status, celery_status]):
        overall_status = "DOWN"
    
    # Return combined status
    return {
        "status": overall_status,
        "components": {
            "database": db_status,
            "redis": redis_status,
            "minio": minio_status,
            "celery": celery_status
        }
    }


@health_router.get("/db", summary="Database health check")
async def database_health_check():
    """
    Database-specific health check with detailed connection information.
    
    Returns:
        dict: Database health status details
    """
    # Get database status
    db_status = check_database()
    
    # Return status with additional details
    return db_status


@health_router.get("/system", summary="System health metrics")
async def system_health_check():
    """
    System metrics health check with memory, CPU, and disk usage.
    
    Returns:
        dict: System health metrics
    """
    # Get system metrics
    metrics = get_system_metrics()
    
    # Return metrics with status
    return {
        "status": "UP",
        "metrics": metrics
    }