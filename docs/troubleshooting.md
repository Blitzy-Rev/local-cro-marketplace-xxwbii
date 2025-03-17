# Introduction

This troubleshooting guide provides comprehensive information for diagnosing and resolving common issues with the Molecular Data Management and CRO Integration Platform. The guide is organized by system component and includes diagnostic steps, common issues, and resolution procedures.

The platform consists of multiple containerized components working together to provide molecular data management and CRO integration capabilities. When issues occur, a systematic approach to troubleshooting can help quickly identify and resolve problems.

This guide is intended for system administrators, IT support staff, and technical users responsible for maintaining the platform. It assumes basic familiarity with Docker, command-line interfaces, and the overall architecture of the system.

# General Troubleshooting Approach

When troubleshooting issues with the platform, follow these general steps:

1. **Identify the Problem**:
   - Gather information about the issue (error messages, affected functionality, timing)
   - Determine which component is likely involved
   - Check if the issue is reproducible

2. **Check System Health**:
   - Verify all containers are running: `docker-compose -f infrastructure/docker-compose.yml ps`
   - Check system resource usage: `docker stats`
   - Use health check endpoints: `curl -s http://localhost/api/v1/health/live`

3. **Check Logs**:
   - View container logs: `docker-compose -f infrastructure/docker-compose.yml logs [service_name]`
   - Check application logs in Kibana (if available)
   - Look for error messages or warnings

4. **Isolate the Issue**:
   - Determine if the issue is with a specific component or system-wide
   - Test related functionality to narrow down the problem
   - Check if recent changes might have caused the issue

5. **Apply Resolution**:
   - Follow component-specific troubleshooting steps
   - Apply the appropriate fix
   - Verify the issue is resolved
   - Document the issue and resolution

For critical issues affecting system availability, refer to the Incident Response section for escalation procedures.

## Using the Health Monitoring Script

The platform includes a health monitoring script that can help diagnose system issues:

```bash
./scripts/monitor-health.sh -v
```

This script checks the health of all system components and provides detailed output about their status. Options include:

- `-v`: Verbose output with detailed information
- `-c`: Continuous monitoring mode
- `-i <seconds>`: Set the interval for continuous monitoring

Example output:

```
System Health Check: 2023-06-10 14:30:45

Container Status:
✅ frontend: Running (Up 2 days)
✅ backend: Running (Up 2 days)
✅ postgres: Running (Up 2 days)
✅ redis: Running (Up 2 days)
✅ minio: Running (Up 2 days)
✅ worker: Running (Up 2 days)
✅ nginx: Running (Up 2 days)

API Health:
✅ Liveness: UP (200 OK, 45ms)
✅ Readiness: UP (200 OK, 78ms)
✅ Database: UP (200 OK, 32ms)

Resource Usage:
💻 CPU: 35% (Normal)
🧠 Memory: 4.2GB/8GB (52%, Normal)
💾 Disk: 45GB/100GB (45%, Normal)

Overall System Status: ✅ HEALTHY
```

If any component shows an error, refer to the specific troubleshooting section for that component.

## Accessing Logs

Logs are essential for troubleshooting. Access logs using the following methods:

1. **Docker Container Logs**:
   ```bash
   # View logs for a specific container
   docker-compose -f infrastructure/docker-compose.yml logs [service_name]
   
   # View logs with timestamps
   docker-compose -f infrastructure/docker-compose.yml logs --timestamps [service_name]
   
   # Follow logs in real-time
   docker-compose -f infrastructure/docker-compose.yml logs --follow [service_name]
   
   # View last 100 lines
   docker-compose -f infrastructure/docker-compose.yml logs --tail=100 [service_name]
   ```

2. **Kibana Log Interface** (if deployed):
   - Access Kibana at `http://localhost:5601`
   - Navigate to "Discover" to search and filter logs
   - Use the pre-configured dashboards for common log views

3. **Log Files**:
   - Application logs: `/var/log/molecular-platform/`
   - Container logs: Managed by Docker's logging driver

When examining logs, look for:
- ERROR or WARNING level messages
- Exception stack traces
- Timing information for performance issues
- Authentication or authorization failures
- Database query errors

## Checking System Resources

Resource constraints can cause various issues. Check system resources using:

```bash
# View resource usage for all containers
docker stats

# Check host system resources
top
df -h
free -m
```

Common resource-related issues:

1. **High CPU Usage**:
   - Slow response times
   - Background tasks taking longer than expected
   - Container restarts due to health check failures

2. **Memory Exhaustion**:
   - Container crashes with "Out of Memory" errors
   - System becoming unresponsive
   - Swap thrashing

3. **Disk Space Issues**:
   - Failed writes or database operations
   - Backup failures
   - Log rotation failures

4. **Network Bottlenecks**:
   - Slow file uploads or downloads
   - Timeouts in inter-service communication
   - Connection failures

If resource constraints are identified, consider:
- Scaling up container resources in `docker-compose.yml`
- Cleaning up unused data and containers
- Adding additional storage
- Optimizing resource-intensive operations

# Container and Deployment Issues

Issues with containers and deployment are common starting points for troubleshooting as they affect the entire system.

## Container Not Starting

If containers fail to start or keep restarting, follow these steps:

1. **Check Container Status**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml ps
   ```
   Look for containers in a "restarting" or "exited" state.

2. **View Container Logs**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml logs [service_name]
   ```
   Look for error messages explaining why the container failed to start.

3. **Check Container Exit Code**:
   ```bash
   docker inspect [container_id] | grep "ExitCode"
   ```
   Common exit codes:
   - 0: Normal exit
   - 1: Application error
   - 137: Container was killed (often due to OOM)
   - 139: Segmentation fault

4. **Verify Environment Variables**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml config
   ```
   Ensure all required environment variables are set correctly.

5. **Check Resource Constraints**:
   ```bash
   docker stats
   ```
   Ensure the host has sufficient resources available.

6. **Verify Volume Mounts**:
   ```bash
   docker inspect [container_id] | grep -A 10 "Mounts"
   ```
   Ensure volumes are correctly mounted and accessible.

7. **Check Network Configuration**:
   ```bash
   docker network ls
   docker network inspect [network_name]
   ```
   Verify network configuration and connectivity.

**Common Solutions**:
- Fix environment variable configuration in `.env` files
- Ensure sufficient system resources are available
- Check file permissions on mounted volumes
- Rebuild the container: `docker-compose -f infrastructure/docker-compose.yml build [service_name]`
- Remove and recreate the container: `docker-compose -f infrastructure/docker-compose.yml up -d --force-recreate [service_name]`

## Container Networking Issues

If containers cannot communicate with each other, follow these steps:

1. **Check Network Configuration**:
   ```bash
   docker network ls
   docker network inspect [network_name]
   ```
   Verify that all expected containers are connected to the appropriate networks.

2. **Test Inter-Container Connectivity**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec [service_name] ping [target_service]
   docker-compose -f infrastructure/docker-compose.yml exec [service_name] curl -v [target_service]:[port]
   ```
   This tests if containers can reach each other.

3. **Check DNS Resolution**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec [service_name] nslookup [target_service]
   ```
   Verify that container names are resolving correctly.

4. **Verify Port Mappings**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml port [service_name] [container_port]
   ```
   Check if ports are correctly mapped to the host.

5. **Inspect Network Interfaces**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec [service_name] ip addr
   ```
   Verify network interfaces are configured correctly.

**Common Solutions**:
- Recreate Docker networks: `docker-compose -f infrastructure/docker-compose.yml down && docker-compose -f infrastructure/docker-compose.yml up -d`
- Check for conflicting port mappings in `docker-compose.yml`
- Verify firewall rules are not blocking container communication
- Ensure container names are used for service discovery, not IP addresses
- Check for network namespace isolation issues

## Volume and Storage Issues

If containers have problems with persistent storage, follow these steps:

1. **Check Volume Status**:
   ```bash
   docker volume ls
   docker volume inspect [volume_name]
   ```
   Verify that volumes exist and are correctly configured.

2. **Check Volume Mounts**:
   ```bash
   docker inspect [container_id] | grep -A 10 "Mounts"
   ```
   Ensure volumes are correctly mounted in containers.

3. **Verify File Permissions**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec [service_name] ls -la [directory]
   ```
   Check if file permissions allow the container to read/write.

4. **Check Disk Space**:
   ```bash
   df -h
   ```
   Ensure sufficient disk space is available.

5. **Inspect Volume Data**:
   ```bash
   docker run --rm -v [volume_name]:/data alpine ls -la /data
   ```
   Verify that expected data exists in the volume.

**Common Solutions**:
- Fix file permissions: `docker-compose -f infrastructure/docker-compose.yml exec [service_name] chown -R [user]:[group] [directory]`
- Free up disk space: `docker system prune -f`
- Recreate volumes (caution: data loss): `docker-compose -f infrastructure/docker-compose.yml down -v && docker-compose -f infrastructure/docker-compose.yml up -d`
- Restore from backup if data is corrupted
- Check for disk I/O issues: `iostat -x 1`

## Deployment and Update Issues

If problems occur during deployment or updates, follow these steps:

1. **Check Deployment Logs**:
   ```bash
   ./scripts/deploy.sh -v
   ```
   Look for error messages during the deployment process.

2. **Verify Configuration**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml config
   ```
   Check for configuration errors in Docker Compose files.

3. **Check for Version Conflicts**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec [service_name] [command_to_check_version]
   ```
   Verify that component versions are compatible.

4. **Verify Database Migrations**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend alembic current
   docker-compose -f infrastructure/docker-compose.yml exec backend alembic history
   ```
   Check if database migrations have been applied correctly.

5. **Check for Stale Containers or Images**:
   ```bash
   docker ps -a
   docker images
   ```
   Look for old containers or images that might be causing conflicts.

**Common Solutions**:
- Roll back to previous version: `./scripts/restore.sh -b [backup_directory]`
- Force rebuild of images: `docker-compose -f infrastructure/docker-compose.yml build --no-cache`
- Clean up Docker environment: `docker system prune -a`
- Manually apply database migrations: `docker-compose -f infrastructure/docker-compose.yml exec backend alembic upgrade head`
- Check for changes in environment variables or configuration files

# Backend API Issues

The Backend API is the core of the system, handling business logic, data processing, and integration between components.

## API Not Responding

If the API is not responding to requests, follow these steps:

1. **Check Container Status**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml ps backend
   ```
   Verify the backend container is running.

2. **Check API Logs**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml logs backend
   ```
   Look for error messages or exceptions.

3. **Verify API Health**:
   ```bash
   curl -s http://localhost/api/v1/health/live
   curl -s http://localhost/api/v1/health/ready
   ```
   Check if the API health endpoints are responding.

4. **Check Database Connectivity**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend python -c "from app.db.session import SessionLocal; db = SessionLocal(); print('Database connection successful');"
   ```
   Verify the API can connect to the database.

5. **Check Resource Usage**:
   ```bash
   docker stats backend
   ```
   Ensure the container has sufficient resources.

6. **Verify Network Connectivity**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend ping postgres
   docker-compose -f infrastructure/docker-compose.yml exec backend ping redis
   docker-compose -f infrastructure/docker-compose.yml exec backend ping minio
   ```
   Check connectivity to dependent services.

**Common Solutions**:
- Restart the API container: `docker-compose -f infrastructure/docker-compose.yml restart backend`
- Increase container resources in `docker-compose.yml`
- Check for deadlocks or infinite loops in application code
- Verify environment variables are correctly set
- Check for database connection exhaustion
- Look for blocked threads or event loop issues

## API Errors

If the API is returning errors, follow these steps:

1. **Check API Logs for Exceptions**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml logs backend
   ```
   Look for exception stack traces and error messages.

2. **Examine Error Responses**:
   ```bash
   curl -v http://localhost/api/v1/[endpoint]
   ```
   Check the HTTP status code and error message in the response.

3. **Verify Request Format**:
   ```bash
   curl -v -X POST -H "Content-Type: application/json" -d '{"key": "value"}' http://localhost/api/v1/[endpoint]
   ```
   Ensure requests are properly formatted.

4. **Check Authentication**:
   ```bash
   curl -v -H "Authorization: Bearer [token]" http://localhost/api/v1/[endpoint]
   ```
   Verify authentication tokens are valid and properly formatted.

5. **Test with Simplified Requests**:
   ```bash
   curl -v http://localhost/api/v1/health/live
   ```
   Start with simple endpoints to isolate the issue.

**Common Error Types and Solutions**:

- **400 Bad Request**:
  - Check request payload format and required fields
  - Verify content type headers
  - Look for validation errors in the response

- **401 Unauthorized**:
  - Verify JWT token is valid and not expired
  - Check authentication configuration
  - Ensure user credentials are correct

- **403 Forbidden**:
  - Check user permissions for the requested resource
  - Verify role-based access control settings
  - Ensure the user has the necessary role

- **404 Not Found**:
  - Verify the endpoint URL is correct
  - Check if the requested resource exists
  - Ensure API routes are correctly configured

- **422 Unprocessable Entity**:
  - Check for business rule violations
  - Verify data consistency requirements
  - Look for detailed error messages in the response

- **500 Internal Server Error**:
  - Check API logs for exception details
  - Look for database errors or connection issues
  - Verify external service dependencies
  - Check for code bugs or edge cases

## Performance Issues

If the API is responding slowly, follow these steps:

1. **Check Response Times**:
   ```bash
   time curl -s http://localhost/api/v1/[endpoint] > /dev/null
   ```
   Measure response times for different endpoints.

2. **Monitor Resource Usage**:
   ```bash
   docker stats backend
   ```
   Check CPU, memory, and I/O usage during requests.

3. **Examine Database Performance**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "SELECT * FROM pg_stat_activity;"
   ```
   Look for slow queries or connection issues.

4. **Check for N+1 Query Problems**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml logs backend | grep "SELECT"
   ```
   Look for repeated similar queries that could be optimized.

5. **Verify Caching**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec redis redis-cli --stat
   ```
   Check if Redis cache is being utilized effectively.

6. **Profile API Endpoints**:
   Enable profiling for specific endpoints and analyze the results.

**Common Solutions**:
- Optimize database queries with proper indexing
- Implement or fix caching for frequently accessed data
- Increase container resources in `docker-compose.yml`
- Optimize expensive operations in application code
- Use batch processing for large datasets
- Implement pagination for large result sets
- Check for blocking operations in asynchronous code

## Molecular Processing Issues

If problems occur with molecular data processing, follow these steps:

1. **Check RDKit Functionality**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend python -c "from rdkit import Chem; mol = Chem.MolFromSmiles('CCO'); print(mol is not None)"
   ```
   Verify that RDKit is working correctly.

2. **Test SMILES Validation**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend python -c "from app.molecular.validator import validate_smiles; print(validate_smiles('CCO')); print(validate_smiles('invalid'))"
   ```
   Check if SMILES validation is working.

3. **Verify Property Calculation**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend python -c "from app.molecular.property_calculator import calculate_properties; from rdkit import Chem; mol = Chem.MolFromSmiles('CCO'); print(calculate_properties(mol))"
   ```
   Test property calculation functionality.

4. **Check for Invalid Molecular Data**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend python -c "from app.db.session import SessionLocal; from app.models.molecule import Molecule; db = SessionLocal(); invalid = db.query(Molecule).filter(Molecule.smiles.like('%invalid%')).all(); print(f'Found {len(invalid)} potentially invalid molecules')"
   ```
   Look for potentially invalid molecular data in the database.

**Common Issues and Solutions**:

- **Invalid SMILES Strings**:
  - Validate and clean SMILES data during import
  - Implement better error handling for invalid structures
  - Use canonical SMILES to avoid representation issues

- **Property Calculation Failures**:
  - Check for molecules that cause calculation errors
  - Implement fallbacks for failed calculations
  - Update RDKit version if bugs are encountered

- **Performance Issues with Large Datasets**:
  - Process molecules in batches
  - Implement caching for computed properties
  - Use background tasks for intensive calculations

- **Memory Issues**:
  - Limit the number of molecules processed simultaneously
  - Optimize memory usage in RDKit operations
  - Increase container memory limits

# Database Issues

The PostgreSQL database stores all structured data for the platform. Database issues can affect multiple system components.

## Connection Issues

If services cannot connect to the database, follow these steps:

1. **Check Database Container Status**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml ps postgres
   ```
   Verify the database container is running.

2. **Check Database Logs**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml logs postgres
   ```
   Look for error messages or connection issues.

3. **Verify Database Connectivity**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres pg_isready
   ```
   Check if the database is accepting connections.

4. **Test Connection from Services**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend python -c "from app.db.session import SessionLocal; db = SessionLocal(); print('Database connection successful');"
   ```
   Verify services can connect to the database.

5. **Check Connection Configuration**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend env | grep POSTGRES
   ```
   Verify database connection parameters.

6. **Check Connection Limits**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -c "SHOW max_connections;"
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "SELECT count(*) FROM pg_stat_activity;"
   ```
   Check if connection limits have been reached.

**Common Solutions**:
- Restart the database container: `docker-compose -f infrastructure/docker-compose.yml restart postgres`
- Verify database credentials in environment variables
- Check network connectivity between services and database
- Increase max_connections if limits are reached
- Check for connection leaks in application code
- Verify database volume permissions

## Data Integrity Issues

If data integrity problems occur, follow these steps:

1. **Check for Constraint Violations**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "SELECT * FROM pg_constraint WHERE contype = 'f' AND conrelid = 'table_name'::regclass;"
   ```
   Identify foreign key constraints for the affected table.

2. **Look for Orphaned Records**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "SELECT a.* FROM table_a a LEFT JOIN table_b b ON a.b_id = b.id WHERE b.id IS NULL;"
   ```
   Find records with missing related records.

3. **Check for Duplicate Data**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "SELECT column, COUNT(*) FROM table GROUP BY column HAVING COUNT(*) > 1;"
   ```
   Identify duplicate values in columns that should be unique.

4. **Verify Data Consistency**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "SELECT * FROM table WHERE condition_that_should_not_exist;"
   ```
   Check for data that violates business rules.

5. **Examine Transaction Logs**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "SELECT * FROM pg_stat_activity WHERE state = 'active';"
   ```
   Look for long-running transactions that might be causing issues.

**Common Solutions**:
- Fix data inconsistencies with targeted SQL updates
- Restore from backup if data corruption is severe
- Implement additional data validation in application code
- Add database constraints to prevent future issues
- Check for race conditions in concurrent operations
- Verify transaction boundaries in application code

## Performance Issues

If database performance is poor, follow these steps:

1. **Identify Slow Queries**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "SELECT query, calls, total_time, mean_time FROM pg_stat_statements ORDER BY mean_time DESC LIMIT 10;"
   ```
   Find the slowest queries in the system.

2. **Check Index Usage**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "SELECT relname, seq_scan, idx_scan FROM pg_stat_user_tables ORDER BY seq_scan DESC;"
   ```
   Identify tables with high sequential scan counts (missing indexes).

3. **Analyze Query Plans**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "EXPLAIN ANALYZE [slow_query];"
   ```
   Examine execution plans for slow queries.

4. **Check Table Statistics**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "SELECT relname, n_live_tup FROM pg_stat_user_tables ORDER BY n_live_tup DESC;"
   ```
   Identify large tables that might need optimization.

5. **Monitor Database Resource Usage**:
   ```bash
   docker stats postgres
   ```
   Check CPU, memory, and I/O usage.

6. **Verify Autovacuum Settings**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "SELECT * FROM pg_settings WHERE name LIKE 'autovacuum%';"
   ```
   Check autovacuum configuration.

**Common Solutions**:
- Add indexes for frequently queried columns
- Optimize complex queries with better join strategies
- Update table statistics: `ANALYZE table_name;`
- Increase shared_buffers for better caching
- Configure autovacuum for more aggressive cleanup
- Implement query timeouts for long-running queries
- Consider table partitioning for very large tables
- Use connection pooling to reduce connection overhead

## Migration Issues

If database migrations fail, follow these steps:

1. **Check Migration Status**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend alembic current
   ```
   Identify the current migration version.

2. **View Migration History**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend alembic history
   ```
   See the full migration history.

3. **Check Migration Logs**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml logs backend | grep -i migration
   ```
   Look for migration-related error messages.

4. **Verify Database Schema**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "\\dt"
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "\\d+ table_name"
   ```
   Examine the current database schema.

5. **Check for Conflicting Changes**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend alembic show [revision]
   ```
   Review specific migration revisions for conflicts.

**Common Solutions**:
- Fix migration scripts to handle edge cases
- Manually apply migrations: `docker-compose -f infrastructure/docker-compose.yml exec backend alembic upgrade head`
- Restore from backup before failed migration
- Create a new migration to fix schema issues
- Manually fix database schema with SQL commands
- Check for data that violates new constraints

# File Storage Issues

The MinIO object storage service manages files such as CSV uploads, molecular images, and experimental results.

## Upload Failures

If file uploads are failing, follow these steps:

1. **Check MinIO Container Status**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml ps minio
   ```
   Verify the MinIO container is running.

2. **Check MinIO Logs**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml logs minio
   ```
   Look for error messages related to uploads.

3. **Verify MinIO Connectivity**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend curl -s http://minio:9000/minio/health/live
   ```
   Check if MinIO is responding to health checks.

4. **Test MinIO Client**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend mc config host add myminio http://minio:9000 ${MINIO_ROOT_USER} ${MINIO_ROOT_PASSWORD}
   docker-compose -f infrastructure/docker-compose.yml exec backend mc ls myminio
   ```
   Verify MinIO client can connect and list buckets.

5. **Check Bucket Existence and Permissions**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend mc ls myminio/bucket_name
   ```
   Verify the target bucket exists and is accessible.

6. **Check Disk Space**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec minio df -h
   ```
   Ensure sufficient disk space is available.

**Common Solutions**:
- Create missing buckets: `docker-compose -f infrastructure/docker-compose.yml exec backend mc mb myminio/bucket_name`
- Fix bucket permissions: `docker-compose -f infrastructure/docker-compose.yml exec backend mc policy set download myminio/bucket_name`
- Restart MinIO: `docker-compose -f infrastructure/docker-compose.yml restart minio`
- Verify MinIO credentials in environment variables
- Check file size limits in application code
- Ensure proper content type is set for uploads
- Free up disk space if storage is full

## Download Failures

If file downloads are failing, follow these steps:

1. **Verify File Existence**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend mc ls myminio/bucket_name/path/to/file
   ```
   Check if the file exists in MinIO.

2. **Test File Access**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend mc cat myminio/bucket_name/path/to/file > /dev/null
   ```
   Verify the file can be accessed.

3. **Check File Metadata**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend mc stat myminio/bucket_name/path/to/file
   ```
   Examine file metadata for issues.

4. **Verify Presigned URL Generation**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend mc share download myminio/bucket_name/path/to/file
   ```
   Test presigned URL generation.

5. **Check Bucket Policies**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend mc policy get myminio/bucket_name
   ```
   Verify bucket access policies.

**Common Solutions**:
- Update bucket policies to allow downloads
- Fix presigned URL generation in application code
- Check for file path encoding issues
- Verify file permissions in MinIO
- Ensure correct content type is set for downloads
- Check for network connectivity issues between services
- Verify that file references in the database match actual files in storage

## Storage Capacity Issues

If storage capacity problems occur, follow these steps:

1. **Check Disk Usage**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec minio df -h
   ```
   Verify available disk space.

2. **Analyze Bucket Usage**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend mc du myminio
   docker-compose -f infrastructure/docker-compose.yml exec backend mc du myminio/bucket_name
   ```
   Identify which buckets are using the most space.

3. **List Large Files**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend mc find myminio/bucket_name --larger-than 10MB
   ```
   Find large files that might be consuming space.

4. **Check for Temporary Files**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend mc ls myminio/temp-files
   ```
   Look for temporary files that should be cleaned up.

5. **Verify Lifecycle Policies**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend mc ilm list myminio/bucket_name
   ```
   Check if lifecycle policies are configured.

**Common Solutions**:
- Implement or fix cleanup jobs for temporary files
- Configure lifecycle policies for automatic cleanup
- Add more storage capacity to the MinIO volume
- Archive old data that is rarely accessed
- Implement file size limits for uploads
- Optimize file formats to reduce storage requirements
- Use compression for large files

## CSV Processing Issues

If CSV file processing is failing, follow these steps:

1. **Check CSV File Existence**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend mc ls myminio/csv-uploads/path/to/file.csv
   ```
   Verify the CSV file was uploaded successfully.

2. **Examine CSV File Content**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend mc cat myminio/csv-uploads/path/to/file.csv | head
   ```
   Check the format and content of the CSV file.

3. **Verify CSV Processing Logs**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml logs backend | grep -i csv
   docker-compose -f infrastructure/docker-compose.yml logs worker | grep -i csv
   ```
   Look for error messages related to CSV processing.

4. **Check Background Task Status**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec redis redis-cli -n 1 keys "celery-task-meta-*"
   docker-compose -f infrastructure/docker-compose.yml exec redis redis-cli -n 1 get "celery-task-meta-[task_id]"
   ```
   Verify the status of CSV processing tasks.

5. **Test CSV Parsing**:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec backend python -c "import pandas as pd; df = pd.read_csv('/tmp/test.csv'); print(df.head())"