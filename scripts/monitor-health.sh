#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Script directory and project root
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(dirname "$SCRIPT_DIR")
INFRASTRUCTURE_DIR="$PROJECT_ROOT/infrastructure"
LOG_DIR="$PROJECT_ROOT/logs"
LOG_FILE="$LOG_DIR/health_$(date +%Y%m%d).log"

# Default configuration
API_HOST="localhost"
API_PORT="80"
API_PROTOCOL="http"
HEALTH_ENDPOINT="/api/health"
TIMEOUT="5"
CHECK_INTERVAL="60"
ALERT_THRESHOLD="3"
ALERT_EMAIL=""
VERBOSE="false"
SINGLE_CHECK="false"
CONTINUOUS="false"
SERVICES="backend postgres redis minio worker"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log() {
  local level="$1"
  local message="$2"
  local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
  
  # Create log directory if it doesn't exist
  mkdir -p "$LOG_DIR"
  
  # Format the message
  local formatted_message="[$timestamp] [$level] $message"
  
  # Output to console with color based on level
  case "$level" in
    "INFO")
      echo -e "${GREEN}$formatted_message${NC}"
      ;;
    "WARNING")
      echo -e "${YELLOW}$formatted_message${NC}"
      ;;
    "ERROR")
      echo -e "${RED}$formatted_message${NC}"
      ;;
    "DEBUG")
      if [[ "$VERBOSE" == "true" ]]; then
        echo -e "${BLUE}$formatted_message${NC}"
      fi
      ;;
    *)
      echo "$formatted_message"
      ;;
  esac
  
  # Append to log file
  echo "$formatted_message" >> "$LOG_FILE"
}

# Print a banner
print_banner() {
  echo -e "${BLUE}"
  echo "┌─────────────────────────────────────────────────────────────┐"
  echo "│ Molecular Data Management and CRO Integration Platform       │"
  echo "│ Health Monitoring System                                    │"
  echo "└─────────────────────────────────────────────────────────────┘"
  echo -e "${NC}"
  echo "Starting health check at $(date)"
  echo "---------------------------------------------------------------"
}

# Check if required dependencies are installed
check_dependencies() {
  local missing_deps=0
  
  # Check for curl
  if ! command -v curl &> /dev/null; then
    log "ERROR" "curl is not installed. Please install curl to use this script."
    missing_deps=1
  fi
  
  # Check for jq
  if ! command -v jq &> /dev/null; then
    log "ERROR" "jq is not installed. Please install jq to use this script."
    missing_deps=1
  fi
  
  # Check for docker
  if ! command -v docker &> /dev/null; then
    log "ERROR" "docker is not installed. Please install docker to use this script."
    missing_deps=1
  fi
  
  # Check for docker-compose
  if ! command -v docker-compose &> /dev/null; then
    log "ERROR" "docker-compose is not installed. Please install docker-compose to use this script."
    missing_deps=1
  fi
  
  if [[ $missing_deps -eq 0 ]]; then
    log "INFO" "All dependencies are satisfied."
    return 0
  else
    return 1
  fi
}

# Check basic liveness of the backend service
check_liveness() {
  log "INFO" "Checking liveness..."
  
  # Make HTTP request to liveness endpoint
  local response
  local status_code
  
  response=$(curl -s -o /dev/null -w "%{http_code}" -m "$TIMEOUT" "$API_PROTOCOL://$API_HOST:$API_PORT$HEALTH_ENDPOINT/live")
  status_code=$?
  
  if [[ $status_code -ne 0 ]]; then
    log "ERROR" "Liveness check failed: Connection error"
    return 1
  elif [[ $response -ne 200 ]]; then
    log "ERROR" "Liveness check failed: HTTP $response"
    return 1
  else
    log "INFO" "Liveness check passed"
    return 0
  fi
}

# Check readiness of the backend service and its dependencies
check_readiness() {
  log "INFO" "Checking readiness..."
  
  # Make HTTP request to readiness endpoint
  local response
  local status_code
  
  response=$(curl -s -m "$TIMEOUT" "$API_PROTOCOL://$API_HOST:$API_PORT$HEALTH_ENDPOINT/ready")
  status_code=$?
  
  if [[ $status_code -ne 0 ]]; then
    log "ERROR" "Readiness check failed: Connection error"
    return 1
  fi
  
  # Parse JSON response
  local status
  status=$(echo "$response" | jq -r '.status')
  
  if [[ "$status" != "UP" ]]; then
    log "ERROR" "Readiness check failed: System not ready"
    
    # Get detailed component statuses
    local db_status=$(echo "$response" | jq -r '.details.database.status')
    local redis_status=$(echo "$response" | jq -r '.details.redis.status')
    
    [[ "$db_status" != "UP" ]] && log "ERROR" "Database status: DOWN"
    [[ "$redis_status" != "UP" ]] && log "ERROR" "Redis status: DOWN"
    
    return 1
  else
    log "INFO" "Readiness check passed: System is ready"
    return 0
  fi
}

# Perform a deep health check of all components
check_deep_health() {
  log "INFO" "Performing deep health check..."
  
  # Make HTTP request to deep health endpoint
  local response
  local status_code
  
  response=$(curl -s -m "$TIMEOUT" "$API_PROTOCOL://$API_HOST:$API_PORT$HEALTH_ENDPOINT/deep")
  status_code=$?
  
  if [[ $status_code -ne 0 ]]; then
    log "ERROR" "Deep health check failed: Connection error"
    return 1
  fi
  
  # Parse JSON response
  local status
  status=$(echo "$response" | jq -r '.status')
  
  if [[ "$status" != "UP" ]]; then
    log "ERROR" "Deep health check failed: System not healthy"
    
    # Get detailed component statuses
    local components=$(echo "$response" | jq -r '.components | keys[]')
    
    for component in $components; do
      local component_status=$(echo "$response" | jq -r ".components.$component.status")
      
      if [[ "$component_status" != "UP" ]]; then
        log "ERROR" "$component status: DOWN"
        local error=$(echo "$response" | jq -r ".components.$component.error // \"Unknown error\"")
        log "ERROR" "$component error: $error"
      else
        local latency=$(echo "$response" | jq -r ".components.$component.latency_ms // \"N/A\"")
        log "INFO" "$component status: UP (latency: $latency ms)"
      fi
    done
    
    return 1
  else
    log "INFO" "Deep health check passed: All components are healthy"
    
    # Show detailed latency information
    local components=$(echo "$response" | jq -r '.components | keys[]')
    
    for component in $components; do
      local latency=$(echo "$response" | jq -r ".components.$component.latency_ms // \"N/A\"")
      log "DEBUG" "$component latency: $latency ms"
    done
    
    return 0
  fi
}

# Check system resource metrics
check_system_health() {
  log "INFO" "Checking system resource metrics..."
  
  # Make HTTP request to system health endpoint
  local response
  local status_code
  
  response=$(curl -s -m "$TIMEOUT" "$API_PROTOCOL://$API_HOST:$API_PORT$HEALTH_ENDPOINT/system")
  status_code=$?
  
  if [[ $status_code -ne 0 ]]; then
    log "ERROR" "System health check failed: Connection error"
    return 1
  fi
  
  # Parse JSON response
  local status
  status=$(echo "$response" | jq -r '.status')
  
  if [[ "$status" != "UP" ]]; then
    log "ERROR" "System health check failed"
    return 1
  else
    log "INFO" "System health check passed"
    
    # Get system metrics
    local cpu_percent=$(echo "$response" | jq -r '.metrics.cpu.percent')
    local memory_percent=$(echo "$response" | jq -r '.metrics.memory.percent')
    local disk_percent=$(echo "$response" | jq -r '.metrics.disk.percent')
    
    log "INFO" "CPU usage: $cpu_percent%"
    log "INFO" "Memory usage: $memory_percent%"
    log "INFO" "Disk usage: $disk_percent%"
    
    # Check for high resource usage warnings
    if [[ $cpu_percent -gt 80 ]]; then
      log "WARNING" "High CPU usage detected"
    fi
    
    if [[ $memory_percent -gt 80 ]]; then
      log "WARNING" "High memory usage detected"
    fi
    
    if [[ $disk_percent -gt 80 ]]; then
      log "WARNING" "High disk usage detected"
    fi
    
    return 0
  fi
}

# Check the status of Docker containers
check_container_status() {
  log "INFO" "Checking container status..."
  
  # Change directory to infrastructure directory where docker-compose.yml is located
  cd "$INFRASTRUCTURE_DIR" || {
    log "ERROR" "Failed to change to infrastructure directory: $INFRASTRUCTURE_DIR"
    return 1
  }
  
  # Get container status
  local container_status
  container_status=$(docker-compose ps)
  
  # Check container state for each service
  local failed_services=0
  
  for service in $SERVICES; do
    if ! echo "$container_status" | grep -q "$service.*Up"; then
      log "ERROR" "$service container is not running"
      failed_services=$((failed_services + 1))
    else
      log "DEBUG" "$service container is running"
    fi
  done
  
  if [[ $failed_services -eq 0 ]]; then
    log "INFO" "All containers are running"
    return 0
  else
    log "ERROR" "$failed_services containers are not running"
    return 1
  fi
}

# Check database connectivity directly
check_database_connection() {
  log "INFO" "Checking direct database connectivity..."
  
  # Change directory to infrastructure directory
  cd "$INFRASTRUCTURE_DIR" || {
    log "ERROR" "Failed to change to infrastructure directory: $INFRASTRUCTURE_DIR"
    return 1
  }
  
  # Use docker-compose exec to run pg_isready
  if docker-compose exec -T postgres pg_isready > /dev/null 2>&1; then
    log "INFO" "Database connection successful"
    return 0
  else
    log "ERROR" "Database connection failed"
    return 1
  fi
}

# Check Redis connectivity directly
check_redis_connection() {
  log "INFO" "Checking direct Redis connectivity..."
  
  # Change directory to infrastructure directory
  cd "$INFRASTRUCTURE_DIR" || {
    log "ERROR" "Failed to change to infrastructure directory: $INFRASTRUCTURE_DIR"
    return 1
  }
  
  # Use docker-compose exec to run redis-cli ping
  local redis_ping
  redis_ping=$(docker-compose exec -T redis redis-cli ping 2>&1)
  
  if [[ "$redis_ping" == "PONG" ]]; then
    log "INFO" "Redis connection successful"
    return 0
  else
    log "ERROR" "Redis connection failed: $redis_ping"
    return 1
  fi
}

# Check MinIO connectivity directly
check_minio_connection() {
  log "INFO" "Checking direct MinIO connectivity..."
  
  # Change directory to infrastructure directory
  cd "$INFRASTRUCTURE_DIR" || {
    log "ERROR" "Failed to change to infrastructure directory: $INFRASTRUCTURE_DIR"
    return 1
  }
  
  # Check MinIO health endpoint
  local minio_health
  minio_health=$(curl -s -m "$TIMEOUT" http://localhost:9000/minio/health/live)
  
  if [[ "$minio_health" == *"Healthy"* ]]; then
    log "INFO" "MinIO connection successful"
    return 0
  else
    log "ERROR" "MinIO connection failed"
    return 1
  fi
}

# Check Celery worker status
check_worker_status() {
  log "INFO" "Checking Celery worker status..."
  
  # Change directory to infrastructure directory
  cd "$INFRASTRUCTURE_DIR" || {
    log "ERROR" "Failed to change to infrastructure directory: $INFRASTRUCTURE_DIR"
    return 1
  }
  
  # Use docker-compose exec to run celery inspect ping
  local worker_ping
  worker_ping=$(docker-compose exec -T worker celery -A app.worker.celery_app inspect ping 2>&1)
  
  if [[ "$worker_ping" == *"OK"* ]]; then
    log "INFO" "Celery worker is active"
    return 0
  else
    log "ERROR" "Celery worker is not responding: $worker_ping"
    return 1
  fi
}

# Run all health checks and return overall status
run_all_checks() {
  local failed_checks=0
  
  # Run all health checks
  check_liveness || failed_checks=$((failed_checks + 1))
  check_readiness || failed_checks=$((failed_checks + 1))
  check_deep_health || failed_checks=$((failed_checks + 1))
  check_system_health || failed_checks=$((failed_checks + 1))
  check_container_status || failed_checks=$((failed_checks + 1))
  check_database_connection || failed_checks=$((failed_checks + 1))
  check_redis_connection || failed_checks=$((failed_checks + 1))
  check_minio_connection || failed_checks=$((failed_checks + 1))
  check_worker_status || failed_checks=$((failed_checks + 1))
  
  # Determine overall health status
  if [[ $failed_checks -eq 0 ]]; then
    log "INFO" "Overall health status: HEALTHY"
    return 0
  else
    local total_checks=9
    local success_rate=$(( (total_checks - failed_checks) * 100 / total_checks ))
    log "ERROR" "Overall health status: UNHEALTHY ($success_rate% checks passed, $failed_checks failures)"
    return 1
  fi
}

# Run a specific health check by name
run_single_check() {
  local check_name="$1"
  
  case "$check_name" in
    "liveness")
      check_liveness
      ;;
    "readiness")
      check_readiness
      ;;
    "deep")
      check_deep_health
      ;;
    "system")
      check_system_health
      ;;
    "containers")
      check_container_status
      ;;
    "database")
      check_database_connection
      ;;
    "redis")
      check_redis_connection
      ;;
    "minio")
      check_minio_connection
      ;;
    "worker")
      check_worker_status
      ;;
    *)
      log "ERROR" "Unknown check name: $check_name"
      log "INFO" "Available checks: liveness, readiness, deep, system, containers, database, redis, minio, worker"
      return 1
      ;;
  esac
  
  return $?
}

# Send an alert when health checks fail repeatedly
send_alert() {
  local message="$1"
  
  # Check if email is configured
  if [[ -n "$ALERT_EMAIL" ]]; then
    log "INFO" "Sending alert email to $ALERT_EMAIL"
    
    # Use mail command if available
    if command -v mail &> /dev/null; then
      echo "$message" | mail -s "ALERT: Molecular Platform Health Check Failure" "$ALERT_EMAIL"
      log "INFO" "Alert email sent"
      return 0
    else
      log "WARNING" "mail command not found, cannot send email alert"
      return 1
    fi
  else
    log "WARNING" "No alert email configured, skipping alert"
    return 1
  fi
}

# Run health checks continuously at specified interval
continuous_monitoring() {
  log "INFO" "Starting continuous monitoring with interval of $CHECK_INTERVAL seconds"
  
  # Initialize failure counter
  local consecutive_failures=0
  
  # Continuous monitoring loop
  while true; do
    # Run all health checks
    if run_all_checks; then
      # Reset failure counter on success
      consecutive_failures=0
      log "INFO" "Health check passed, no alerts needed"
    else
      # Increment failure counter
      consecutive_failures=$((consecutive_failures + 1))
      log "WARNING" "Health check failed, consecutive failures: $consecutive_failures"
      
      # Send alert if threshold is reached
      if [[ $consecutive_failures -ge $ALERT_THRESHOLD ]]; then
        log "ERROR" "Alert threshold reached ($consecutive_failures failures)"
        local alert_message="The Molecular Platform health check has failed $consecutive_failures consecutive times.\n\nPlease check the system immediately.\n\nTimestamp: $(date)\n\nLog file: $LOG_FILE"
        send_alert "$alert_message"
        # Reset counter after alert is sent
        consecutive_failures=0
      fi
    fi
    
    # Wait for next check interval
    log "DEBUG" "Sleeping for $CHECK_INTERVAL seconds until next check"
    sleep "$CHECK_INTERVAL"
  done
}

# Parse command line arguments
parse_arguments() {
  local OPTIND
  
  # Process options
  while getopts ":h:p:e:t:i:a:c:vsl:C" opt; do
    case "$opt" in
      h)
        API_HOST="$OPTARG"
        ;;
      p)
        API_PORT="$OPTARG"
        ;;
      e)
        HEALTH_ENDPOINT="$OPTARG"
        ;;
      t)
        TIMEOUT="$OPTARG"
        ;;
      i)
        CHECK_INTERVAL="$OPTARG"
        ;;
      a)
        ALERT_THRESHOLD="$OPTARG"
        ;;
      c)
        ALERT_EMAIL="$OPTARG"
        ;;
      v)
        VERBOSE="true"
        ;;
      s)
        SINGLE_CHECK="true"
        ;;
      l)
        # Single check with specified name
        SINGLE_CHECK="true"
        SINGLE_CHECK_NAME="$OPTARG"
        ;;
      C)
        CONTINUOUS="true"
        ;;
      \?)
        echo "Invalid option: -$OPTARG" >&2
        display_usage
        exit 1
        ;;
      :)
        echo "Option -$OPTARG requires an argument." >&2
        display_usage
        exit 1
        ;;
    esac
  done
  
  # Shift to process remaining arguments
  shift $((OPTIND-1))
  
  # If we're in single check mode and no check name was provided, use the first argument if available
  if [[ "$SINGLE_CHECK" == "true" && -z "$SINGLE_CHECK_NAME" && -n "$1" ]]; then
    SINGLE_CHECK_NAME="$1"
  fi
}

# Display usage information
display_usage() {
  echo "Usage: $0 [OPTIONS] [CHECK_NAME]"
  echo
  echo "Monitor the health of the Molecular Data Management and CRO Integration Platform."
  echo
  echo "Options:"
  echo "  -h HOST         API host (default: localhost)"
  echo "  -p PORT         API port (default: 80)"
  echo "  -e ENDPOINT     Base health endpoint (default: /api/health)"
  echo "  -t SECONDS      Request timeout in seconds (default: 5)"
  echo "  -i SECONDS      Check interval for continuous monitoring (default: 60)"
  echo "  -a COUNT        Alert threshold for consecutive failures (default: 3)"
  echo "  -c EMAIL        Email address for alerts"
  echo "  -v              Verbose output (include debug messages)"
  echo "  -s              Single check mode (run once and exit)"
  echo "  -l CHECK_NAME   Run a specific check (liveness, readiness, deep, system, containers, database, redis, minio, worker)"
  echo "  -C              Continuous monitoring mode"
  echo
  echo "CHECK_NAME (optional):"
  echo "  When used with -l, specifies which check to run"
  echo
  echo "Examples:"
  echo "  $0 -s                       # Run all checks once"
  echo "  $0 -l deep                  # Run only the deep health check"
  echo "  $0 -C -i 300 -a 5 -c admin@example.com  # Continuous monitoring with alerts"
  echo
}

# Main function that orchestrates the health monitoring process
main() {
  # Parse command line arguments
  parse_arguments "$@"
  
  # Display banner unless in single check mode
  if [[ "$SINGLE_CHECK" != "true" ]]; then
    print_banner
  fi
  
  # Check dependencies
  check_dependencies || exit 1
  
  # Run in the appropriate mode
  if [[ "$SINGLE_CHECK" == "true" && -n "$SINGLE_CHECK_NAME" ]]; then
    # Run a specific check
    log "INFO" "Running single check: $SINGLE_CHECK_NAME"
    run_single_check "$SINGLE_CHECK_NAME"
    exit $?
  elif [[ "$CONTINUOUS" == "true" ]]; then
    # Run continuous monitoring
    continuous_monitoring
  else
    # Run all checks once
    run_all_checks
    exit $?
  fi
}

# Execute main function with all arguments
main "$@"