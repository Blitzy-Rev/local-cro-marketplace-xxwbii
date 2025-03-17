# First stage of the multi-stage build process that installs dependencies and prepares the application
FROM python:3.10-slim AS builder
# Set Python environment variables to optimize container behavior
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
# Set the working directory for subsequent instructions
WORKDIR /app
# Install system dependencies required for building Python packages and RDKit
RUN apt-get update && apt-get install -y --no-install-recommends build-essential libpq-dev curl git gcc g++ cmake libboost-all-dev libeigen3-dev libfreetype6-dev libpng-dev && apt-get clean && rm -rf /var/lib/apt/lists/*
# Copy dependency definition file
COPY requirements.txt /app/requirements.txt
# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt
# Copy the application code
COPY . /app/

# Start the final stage with a clean base image
FROM python:3.10-slim
# Install minimal runtime dependencies including RDKit requirements
RUN apt-get update && apt-get install -y --no-install-recommends libpq5 libboost-all-dev libeigen3-dev libfreetype6-dev libpng-dev && apt-get clean && rm -rf /var/lib/apt/lists/*
# Copy installed packages from builder stage
COPY --from=builder /usr/local/lib/python3.10/site-packages /usr/local/lib/python3.10/site-packages
# Copy application code from builder stage
COPY --from=builder /app /app
# Copy the entrypoint script
COPY entrypoint.sh /entrypoint.sh
# Make entrypoint script executable
RUN chmod +x /entrypoint.sh
# Set environment variables for the application
ENV PYTHONPATH=/app PORT=8000 APP_MODULE=main:app WORKERS_PER_CORE=1 MAX_WORKERS=4 LOG_LEVEL=info
# Set the working directory
WORKDIR /app
# Expose the port the application will run on
EXPOSE 8000
# Set the entrypoint script
ENTRYPOINT ["/entrypoint.sh"]
# Set the default command to run the application with Gunicorn and Uvicorn workers
CMD ["gunicorn", "-k", "uvicorn.workers.UvicornWorker", "-c", "/app/gunicorn_conf.py", "main:app"]