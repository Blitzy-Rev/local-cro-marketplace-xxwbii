from datetime import datetime, timedelta
from typing import Dict, List, Optional, Any, Tuple

import os
import sys
import platform
import psutil  # psutil version 5.9.0+

from sqlalchemy import select, func, desc  # sqlalchemy version 2.0+
from sqlalchemy.orm import Session  # sqlalchemy.orm version 2.0+

from ..crud.crud_user import user
from ..models.user import User
from ..models.molecule import Molecule
from ..models.library import Library
from ..models.experiment import Experiment
from ..models.submission import Submission
from ..models.result import Result
from ..models.notification import Notification
from .file_storage_service import FileStorageService
from ..core.security import get_password_hash
from ..exceptions import AdminException
from ..logging_config import logger


class AdminService:
    """
    Service for administrative functions including user management,
    system monitoring, and activity logging.
    """
    
    def __init__(self, db: Session, file_storage_service: FileStorageService):
        """
        Initialize the admin service with database session and file storage service.
        
        Args:
            db: Database session
            file_storage_service: File storage service for system monitoring
        """
        self.db = db
        self.file_storage_service = file_storage_service
    
    def get_users(self, role: Optional[str] = None, status: Optional[str] = None, 
                  skip: int = 0, limit: int = 100) -> List[User]:
        """
        Get list of users with optional filtering by role and status.
        
        Args:
            role: Filter by user role (pharma, cro, admin)
            status: Filter by user status (active, inactive, pending, locked)
            skip: Number of records to skip for pagination
            limit: Maximum number of records to return
            
        Returns:
            List of user objects matching the criteria
        """
        # If role is provided, use role-specific query
        if role:
            users = user.get_multi_by_role(self.db, role=role, skip=skip, limit=limit)
            
            # Further filter by status if provided
            if status and users:
                users = [u for u in users if u.status.name.lower() == status.lower()]
                
            return users
        
        # Get all users and filter by status if needed
        users = user.get_multi(self.db, skip=skip, limit=limit)
        
        if status and users:
            users = [u for u in users if u.status.name.lower() == status.lower()]
            
        return users
    
    def get_user(self, user_id: int) -> User:
        """
        Get a user by ID.
        
        Args:
            user_id: ID of the user to retrieve
            
        Returns:
            User object if found
            
        Raises:
            AdminException: If user not found
        """
        db_user = user.get(self.db, id=user_id)
        if not db_user:
            raise AdminException(f"User with ID {user_id} not found",
                               details={"user_id": user_id})
        return db_user
    
    def create_user(self, user_data: Dict[str, Any]) -> User:
        """
        Create a new user with the provided data.
        
        Args:
            user_data: Dictionary containing user data (email, password, role)
            
        Returns:
            Created user object
            
        Raises:
            AdminException: If user with email already exists
        """
        # Check if user with email already exists
        email = user_data.get("email")
        existing_user = user.get_by_email(self.db, email=email)
        
        if existing_user:
            raise AdminException(f"User with email {email} already exists",
                               details={"email": email})
        
        # Create new user
        new_user = user.create_with_password(self.db, obj_in=user_data)
        
        # Log user creation
        logger.info(f"Created user with email {email} and role {user_data.get('role')}")
        
        return new_user
    
    def update_user(self, user_id: int, user_data: Dict[str, Any]) -> User:
        """
        Update user information.
        
        Args:
            user_id: ID of the user to update
            user_data: Dictionary containing user data to update
            
        Returns:
            Updated user object
            
        Raises:
            AdminException: If user not found or email already exists for another user
        """
        # Get user to update
        db_user = self.get_user(user_id)
        
        # Check if trying to change email to one that already exists
        if "email" in user_data and user_data["email"] != db_user.email:
            existing_user = user.get_by_email(self.db, email=user_data["email"])
            if existing_user and existing_user.id != user_id:
                raise AdminException(
                    f"User with email {user_data['email']} already exists",
                    details={"email": user_data["email"]}
                )
        
        # Update user
        updated_user = user.update_with_password(self.db, db_obj=db_user, obj_in=user_data)
        
        # Log user update
        logger.info(f"Updated user with ID {user_id}")
        
        return updated_user
    
    def delete_user(self, user_id: int) -> bool:
        """
        Delete a user by ID.
        
        Args:
            user_id: ID of the user to delete
            
        Returns:
            True if user was deleted successfully
            
        Raises:
            AdminException: If user not found
        """
        # Get user to delete (will raise exception if not found)
        db_user = self.get_user(user_id)
        
        # Delete user
        user.remove(self.db, id=user_id)
        
        # Log user deletion
        logger.info(f"Deleted user with ID {user_id} and email {db_user.email}")
        
        return True
    
    def reset_user_password(self, user_id: int, new_password: str) -> User:
        """
        Reset a user's password.
        
        Args:
            user_id: ID of the user
            new_password: New password to set
            
        Returns:
            Updated user object
            
        Raises:
            AdminException: If user not found
        """
        # Get user (will raise exception if not found)
        db_user = self.get_user(user_id)
        
        # Update user with new password
        updated_user = user.update_with_password(
            self.db, 
            db_obj=db_user, 
            obj_in={"password": new_password}
        )
        
        # Log password reset
        logger.info(f"Reset password for user with ID {user_id}")
        
        return updated_user
    
    def get_system_stats(self) -> Dict[str, Any]:
        """
        Get system statistics including user counts, data counts, and resource usage.
        
        Returns:
            Dictionary containing system statistics
        """
        # Get user counts
        total_users = self.db.query(User).count()
        pharma_users = user.count_by_role(self.db, role='pharma')
        cro_users = user.count_by_role(self.db, role='cro')
        admin_users = user.count_by_role(self.db, role='admin')
        
        # Get active users count
        stmt = select(func.count(User.id)).where(User.status == 'ACTIVE')
        active_users = self.db.execute(stmt).scalar_one_or_none() or 0
        
        # Get data counts
        total_molecules = self.db.query(Molecule).count()
        total_libraries = self.db.query(Library).count()
        total_experiments = self.db.query(Experiment).count()
        total_submissions = self.db.query(Submission).count()
        total_results = self.db.query(Result).count()
        
        # Get resource usage
        resource_usage = self._get_resource_usage()
        
        # Compile statistics
        stats = {
            "timestamp": datetime.utcnow(),
            "users": {
                "total": total_users,
                "active": active_users,
                "by_role": {
                    "pharma": pharma_users,
                    "cro": cro_users,
                    "admin": admin_users
                }
            },
            "data": {
                "molecules": total_molecules,
                "libraries": total_libraries,
                "experiments": total_experiments,
                "submissions": total_submissions,
                "results": total_results
            },
            "resources": resource_usage
        }
        
        return stats
    
    def get_system_resources(self) -> Dict[str, Any]:
        """
        Get detailed system resource usage information.
        
        Returns:
            Dictionary containing detailed resource information
        """
        # Get basic resource usage
        basic_usage = self._get_resource_usage()
        
        # Get more detailed information
        # CPU details
        cpu_info = {
            "logical_cores": psutil.cpu_count(logical=True),
            "physical_cores": psutil.cpu_count(logical=False),
            "usage_per_core": psutil.cpu_percent(interval=0.1, percpu=True),
            "load_avg": psutil.getloadavg() if hasattr(psutil, 'getloadavg') else None
        }
        
        # Memory details
        virtual_memory = psutil.virtual_memory()
        memory_info = {
            "total": virtual_memory.total,
            "available": virtual_memory.available,
            "used": virtual_memory.used,
            "percent": virtual_memory.percent,
        }
        
        # Disk details
        disk_info = {}
        for partition in psutil.disk_partitions():
            try:
                usage = psutil.disk_usage(partition.mountpoint)
                disk_info[partition.mountpoint] = {
                    "total": usage.total,
                    "used": usage.used,
                    "free": usage.free,
                    "percent": usage.percent,
                    "device": partition.device,
                    "fstype": partition.fstype
                }
            except (PermissionError, FileNotFoundError):
                continue
        
        # Network details
        try:
            network_info = psutil.net_io_counters(pernic=True)
            network_stats = {
                interface: {
                    "bytes_sent": stats.bytes_sent,
                    "bytes_received": stats.bytes_recv,
                    "packets_sent": stats.packets_sent,
                    "packets_received": stats.packets_recv
                }
                for interface, stats in network_info.items()
            }
        except Exception as e:
            logger.error(f"Error collecting network statistics: {str(e)}")
            network_stats = {"error": str(e)}
        
        # Process information
        process_count = len(psutil.pids())
        
        # Database statistics
        try:
            # Query database size
            db_size_query = """
            SELECT pg_database_size(current_database()) as size;
            """
            db_size = self.db.execute(db_size_query).scalar_one()
            
            # Query connection count
            connection_query = """
            SELECT count(*) FROM pg_stat_activity WHERE datname = current_database();
            """
            connection_count = self.db.execute(connection_query).scalar_one()
            
            db_stats = {
                "size_bytes": db_size,
                "connection_count": connection_count
            }
        except Exception as e:
            logger.error(f"Error collecting database statistics: {str(e)}")
            db_stats = {"error": str(e)}
        
        # Storage statistics
        try:
            storage_stats = {}
            for bucket_name in ["csv-uploads", "molecule-images", "experiment-files", "result-files"]:
                try:
                    files = self.file_storage_service.list_files(bucket_name)
                    total_size = sum(file.get("size", 0) for file in files)
                    storage_stats[bucket_name] = {
                        "file_count": len(files),
                        "total_size_bytes": total_size
                    }
                except Exception as bucket_error:
                    storage_stats[bucket_name] = {"error": str(bucket_error)}
        except Exception as e:
            logger.error(f"Error collecting storage statistics: {str(e)}")
            storage_stats = {"error": str(e)}
        
        # Compile all resource information
        resources = {
            "timestamp": datetime.utcnow(),
            "system": {
                "platform": platform.platform(),
                "processor": platform.processor(),
                "python_version": platform.python_version(),
                "hostname": platform.node()
            },
            "cpu": cpu_info,
            "memory": memory_info,
            "disk": disk_info,
            "network": network_stats,
            "processes": {
                "count": process_count
            },
            "database": db_stats,
            "storage": storage_stats,
            **basic_usage
        }
        
        return resources
    
    def check_system_health(self) -> Dict[str, Any]:
        """
        Perform a comprehensive system health check.
        
        Returns:
            Dictionary containing health check results
        """
        health_results = {
            "timestamp": datetime.utcnow(),
            "status": "UP",  # Default overall status
            "components": {}
        }
        
        # Check database connectivity
        try:
            # Simple query to check database connection
            self.db.execute(select(func.count(User.id)))
            health_results["components"]["database"] = {
                "status": "UP",
                "message": "Database connection successful"
            }
        except Exception as e:
            health_results["components"]["database"] = {
                "status": "DOWN",
                "message": f"Database connection failed: {str(e)}"
            }
            health_results["status"] = "DOWN"  # Critical component
        
        # Check file storage connectivity
        try:
            # List files in a bucket to check storage connection
            self.file_storage_service.list_files("csv-uploads")
            health_results["components"]["file_storage"] = {
                "status": "UP",
                "message": "File storage connection successful"
            }
        except Exception as e:
            health_results["components"]["file_storage"] = {
                "status": "DOWN",
                "message": f"File storage connection failed: {str(e)}"
            }
            health_results["status"] = "DOWN"  # Critical component
        
        # Check resource usage
        try:
            resource_usage = self._get_resource_usage()
            
            # CPU check
            cpu_status = "UP"
            cpu_message = "CPU usage is normal"
            
            if resource_usage["cpu_percent"] > 90:
                cpu_status = "DOWN"
                cpu_message = f"Critical CPU usage: {resource_usage['cpu_percent']}%"
                health_results["status"] = "DOWN"  # Critical component
            elif resource_usage["cpu_percent"] > 75:
                cpu_status = "DEGRADED"
                cpu_message = f"High CPU usage: {resource_usage['cpu_percent']}%"
                # Don't downgrade overall status if it's already DOWN
                if health_results["status"] == "UP":
                    health_results["status"] = "DEGRADED"
            
            health_results["components"]["cpu"] = {
                "status": cpu_status,
                "message": cpu_message,
                "metrics": {
                    "usage_percent": resource_usage["cpu_percent"]
                }
            }
            
            # Memory check
            memory_status = "UP"
            memory_message = "Memory usage is normal"
            
            if resource_usage["memory_percent"] > 90:
                memory_status = "DOWN"
                memory_message = f"Critical memory usage: {resource_usage['memory_percent']}%"
                health_results["status"] = "DOWN"  # Critical component
            elif resource_usage["memory_percent"] > 75:
                memory_status = "DEGRADED"
                memory_message = f"High memory usage: {resource_usage['memory_percent']}%"
                # Don't downgrade overall status if it's already DOWN
                if health_results["status"] == "UP":
                    health_results["status"] = "DEGRADED"
            
            health_results["components"]["memory"] = {
                "status": memory_status,
                "message": memory_message,
                "metrics": {
                    "usage_percent": resource_usage["memory_percent"]
                }
            }
            
            # Disk check
            disk_status = "UP"
            disk_message = "Disk usage is normal"
            
            if resource_usage["disk_percent"] > 90:
                disk_status = "DOWN"
                disk_message = f"Critical disk usage: {resource_usage['disk_percent']}%"
                health_results["status"] = "DOWN"  # Critical component
            elif resource_usage["disk_percent"] > 75:
                disk_status = "DEGRADED"
                disk_message = f"High disk usage: {resource_usage['disk_percent']}%"
                # Don't downgrade overall status if it's already DOWN
                if health_results["status"] == "UP":
                    health_results["status"] = "DEGRADED"
            
            health_results["components"]["disk"] = {
                "status": disk_status,
                "message": disk_message,
                "metrics": {
                    "usage_percent": resource_usage["disk_percent"]
                }
            }
            
        except Exception as e:
            health_results["components"]["resources"] = {
                "status": "UNKNOWN",
                "message": f"Error checking resource usage: {str(e)}"
            }
        
        return health_results
    
    def get_activity_logs(self, filter_params: Dict[str, Any], page: int = 1, size: int = 20) -> Dict[str, Any]:
        """
        Get activity logs with filtering and pagination.
        
        Args:
            filter_params: Dictionary of filter parameters (user_id, action, date_range, etc.)
            page: Page number for pagination (1-based)
            size: Page size for pagination
            
        Returns:
            Dictionary containing paginated activity logs
        """
        # For this implementation, we'll use direct SQL queries since the ActivityLog model
        # is not provided in the imports
        
        # Calculate offset for pagination
        offset = (page - 1) * size
        
        # Build base SQL query
        base_query = """
        SELECT * FROM activity_logs 
        WHERE 1=1
        """
        
        # Build count query
        count_query = """
        SELECT COUNT(*) FROM activity_logs 
        WHERE 1=1
        """
        
        # Build query parameters
        params = {}
        
        # Add filters if provided
        if "user_id" in filter_params and filter_params["user_id"]:
            filter_clause = " AND user_id = :user_id"
            base_query += filter_clause
            count_query += filter_clause
            params["user_id"] = filter_params["user_id"]
        
        if "action" in filter_params and filter_params["action"]:
            filter_clause = " AND action = :action"
            base_query += filter_clause
            count_query += filter_clause
            params["action"] = filter_params["action"]
        
        if "resource_type" in filter_params and filter_params["resource_type"]:
            filter_clause = " AND resource_type = :resource_type"
            base_query += filter_clause
            count_query += filter_clause
            params["resource_type"] = filter_params["resource_type"]
        
        if "date_from" in filter_params and filter_params["date_from"]:
            filter_clause = " AND timestamp >= :date_from"
            base_query += filter_clause
            count_query += filter_clause
            params["date_from"] = filter_params["date_from"]
        
        if "date_to" in filter_params and filter_params["date_to"]:
            filter_clause = " AND timestamp <= :date_to"
            base_query += filter_clause
            count_query += filter_clause
            params["date_to"] = filter_params["date_to"]
        
        # Add sorting (newest first by default)
        base_query += " ORDER BY timestamp DESC"
        
        # Add pagination
        base_query += " LIMIT :limit OFFSET :offset"
        params["limit"] = size
        params["offset"] = offset
        
        try:
            # Execute query
            result = self.db.execute(base_query, params)
            logs = [dict(row) for row in result]
            
            # Get total count for pagination
            total = self.db.execute(count_query, params).scalar_one()
            
            # Calculate pagination metadata
            total_pages = (total + size - 1) // size  # Ceiling division
            
            # Prepare response
            response = {
                "items": logs,
                "total": total,
                "page": page,
                "size": size,
                "pages": total_pages,
            }
            
            return response
        except Exception as e:
            logger.error(f"Error retrieving activity logs: {str(e)}")
            # Return empty result on error
            return {
                "items": [],
                "total": 0,
                "page": page,
                "size": size,
                "pages": 0,
                "error": str(e)
            }
    
    def get_system_alerts(self, severity: Optional[str] = None, 
                          alert_type: Optional[str] = None,
                          resolved: Optional[bool] = None) -> Dict[str, Any]:
        """
        Get system alerts with optional filtering.
        
        Args:
            severity: Filter by alert severity (critical, warning, info)
            alert_type: Filter by alert type (system, performance, security, etc.)
            resolved: Filter by resolution status
            
        Returns:
            Dictionary containing system alerts and counts
        """
        # Use the Notification model for system alerts
        stmt = select(Notification).where(Notification.type == 'SYSTEM_ALERT')
        
        # Apply additional filters if provided
        if severity:
            stmt = stmt.where(Notification.data.contains({"severity": severity}))
        
        if alert_type:
            stmt = stmt.where(Notification.data.contains({"alert_type": alert_type}))
        
        if resolved is not None:
            stmt = stmt.where(Notification.read_status == resolved)
        
        # Apply sorting (newest first)
        stmt = stmt.order_by(desc(Notification.created_at))
        
        # Execute query
        result = self.db.execute(stmt).scalars().all()
        
        # Convert to list of dictionaries
        alerts = []
        for notification in result:
            alert = {
                "id": notification.id,
                "message": notification.message,
                "created_at": notification.created_at,
                "resolved": notification.read_status,
                "resolved_at": notification.read_at,
                "severity": notification.data.get("severity") if notification.data else None,
                "alert_type": notification.data.get("alert_type") if notification.data else None,
                "details": notification.data.get("details") if notification.data else None
            }
            alerts.append(alert)
        
        # Get count by severity for summary
        severity_counts = {}
        for alert in alerts:
            severity = alert.get("severity")
            if severity:
                severity_counts[severity] = severity_counts.get(severity, 0) + 1
        
        # Prepare response
        response = {
            "items": alerts,
            "total": len(alerts),
            "counts_by_severity": severity_counts,
        }
        
        return response
    
    def update_system_alert(self, alert_id: int, alert_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update system alert status.
        
        Args:
            alert_id: ID of the alert to update
            alert_data: Dictionary containing update data (resolved, notes, etc.)
            
        Returns:
            Updated system alert
            
        Raises:
            AdminException: If alert not found
        """
        # Get the notification representing the alert
        stmt = select(Notification).where(
            Notification.id == alert_id,
            Notification.type == 'SYSTEM_ALERT'
        )
        notification = self.db.execute(stmt).scalar_one_or_none()
        
        if not notification:
            raise AdminException(f"Alert with ID {alert_id} not found",
                               details={"alert_id": alert_id})
        
        # Update resolution status if provided
        if "resolved" in alert_data:
            notification.read_status = alert_data["resolved"]
            
            # Set resolved_at timestamp if resolving
            if alert_data["resolved"]:
                notification.read_at = datetime.utcnow()
        
        # Update data field with notes if provided
        if "notes" in alert_data and notification.data:
            data = notification.data.copy() if notification.data else {}
            data["notes"] = alert_data["notes"]
            notification.data = data
        
        # Save changes
        self.db.add(notification)
        self.db.commit()
        self.db.refresh(notification)
        
        # Log alert update
        logger.info(f"Updated system alert with ID {alert_id}")
        
        # Convert to response format
        updated_alert = {
            "id": notification.id,
            "message": notification.message,
            "created_at": notification.created_at,
            "resolved": notification.read_status,
            "resolved_at": notification.read_at,
            "severity": notification.data.get("severity") if notification.data else None,
            "alert_type": notification.data.get("alert_type") if notification.data else None,
            "notes": notification.data.get("notes") if notification.data else None,
            "details": notification.data.get("details") if notification.data else None
        }
        
        return updated_alert
    
    def create_activity_log(self, user_id: int, action: str, resource_type: str,
                            resource_id: Optional[int] = None, 
                            details: Optional[Dict[str, Any]] = None,
                            ip_address: Optional[str] = None) -> Dict[str, Any]:
        """
        Create an activity log entry.
        
        Args:
            user_id: ID of the user performing the action
            action: Description of the action performed
            resource_type: Type of resource being acted upon
            resource_id: ID of the resource (optional)
            details: Additional details about the action (optional)
            ip_address: IP address of the user (optional)
            
        Returns:
            Created activity log entry
        """
        # Since we don't have an ActivityLog model in the imports, we'll use direct SQL
        # This approach ensures we don't make assumptions about schema structure
        
        try:
            # Prepare insert data
            insert_query = """
            INSERT INTO activity_logs 
            (user_id, action, resource_type, resource_id, details, ip_address, timestamp)
            VALUES 
            (:user_id, :action, :resource_type, :resource_id, :details, :ip_address, :timestamp)
            RETURNING id
            """
            
            params = {
                "user_id": user_id,
                "action": action,
                "resource_type": resource_type,
                "resource_id": resource_id,
                "details": details,
                "ip_address": ip_address,
                "timestamp": datetime.utcnow()
            }
            
            # Execute insert and get ID
            result = self.db.execute(insert_query, params)
            log_id = result.scalar_one()
            
            # Prepare response
            log_entry = {
                "id": log_id,
                **params
            }
            
            return log_entry
        
        except Exception as e:
            logger.error(f"Error creating activity log: {str(e)}")
            # Return basic info even if insert fails
            return {
                "user_id": user_id,
                "action": action,
                "resource_type": resource_type,
                "resource_id": resource_id,
                "timestamp": datetime.utcnow(),
                "error": str(e)
            }
    
    def create_system_alert(self, title: str, message: str, severity: str, 
                            alert_type: str, details: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Create a system alert.
        
        Args:
            title: Alert title
            message: Alert message
            severity: Alert severity (critical, warning, info)
            alert_type: Alert type (system, performance, security, etc.)
            details: Additional details about the alert (optional)
            
        Returns:
            Created system alert
        """
        # Use the Notification model for system alerts
        alert_data = {
            "severity": severity,
            "alert_type": alert_type,
            "title": title,
            "details": details
        }
        
        # Create notification with system alert type
        notification = Notification(
            user_id=None,  # System alerts don't belong to a specific user
            type='SYSTEM_ALERT',
            message=message,
            data=alert_data,
            read_status=False,
            created_at=datetime.utcnow()
        )
        
        # Save to database
        self.db.add(notification)
        self.db.commit()
        self.db.refresh(notification)
        
        # Convert to response format
        created_alert = {
            "id": notification.id,
            "title": title,
            "message": notification.message,
            "severity": severity,
            "alert_type": alert_type,
            "details": details,
            "created_at": notification.created_at,
            "resolved": False
        }
        
        return created_alert
    
    def _get_resource_usage(self) -> Dict[str, float]:
        """
        Get basic resource usage information (CPU, memory, disk).
        
        Returns:
            Dictionary containing resource usage percentages
        """
        # Get CPU usage
        cpu_percent = psutil.cpu_percent(interval=0.1)
        
        # Get memory usage
        memory = psutil.virtual_memory()
        memory_percent = memory.percent
        
        # Get disk usage for root partition
        disk = psutil.disk_usage('/')
        disk_percent = disk.percent
        
        return {
            "cpu_percent": cpu_percent,
            "memory_percent": memory_percent,
            "disk_percent": disk_percent
        }