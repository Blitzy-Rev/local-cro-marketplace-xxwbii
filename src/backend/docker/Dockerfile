# ----------------------
# Builder stage for dependencies
# ----------------------
FROM python:3.10-slim AS builder

# Set environment variables to optimize Python behavior
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

# Install system dependencies required for building Python packages and RDKit
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libpq-dev \
    curl \
    git \
    gcc \
    g++ \
    cmake \
    libboost-all-dev \
    libeigen3-dev \
    libfreetype6-dev \
    libpng-dev

# Install Python dependencies
COPY ../requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy the application code
COPY ../ /app/

# ----------------------
# Final stage for runtime
# ----------------------
FROM python:3.10-slim

# Install runtime dependencies only
RUN apt-get update && apt-get install -y --no-install-recommends \
    libpq5 \
    libboost-all-dev \
    libeigen3-dev \
    libfreetype6-dev \
    libpng-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy installed packages from builder stage
COPY --from=builder /usr/local/lib/python3.10/site-packages /usr/local/lib/python3.10/site-packages

# Copy application code from builder stage
COPY --from=builder /app /app

# Copy entrypoint and start scripts
COPY entrypoint.sh /entrypoint.sh
COPY start.sh /start.sh

# Make scripts executable
RUN chmod +x /entrypoint.sh /start.sh

# Set environment variables
ENV PYTHONPATH=/app \
    PORT=8000 \
    APP_MODULE=main:app \
    WORKERS_PER_CORE=1 \
    MAX_WORKERS=4 \
    LOG_LEVEL=info

# Set working directory
WORKDIR /app

# Expose the application port
EXPOSE 8000

# Set entrypoint
ENTRYPOINT ["/entrypoint.sh"]

# Set default command
CMD ["/start.sh"]