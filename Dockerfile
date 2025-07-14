 Dockerfile for GenBank to GTF/GFF3 Converter

# Start from a specific Python version (3.11 slim for smaller image size)
FROM python:3.11-slim

# Set environment variables to optimize Python in Docker
# PYTHONDONTWRITEBYTECODE: prevents Python from writing .pyc files to disk
# PYTHONUNBUFFERED: ensures Python output is unbuffered, useful for logging
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# Install system dependencies required by Python packages or the application.
# `build-essential` is often useful for compiling certain Python packages that
# have C extensions (though Biopython often doesn't strictly need it for basic parsing).
# `libgff3-dev` or similar might be needed if you were using C libraries for GFF3,
# but Biopython handles this in Python.
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential && \
    # Clean up the apt cache to keep the image small
    rm -rf /var/lib/apt/lists/*

# Set the working directory in the container
WORKDIR /app

# Copy the requirements file and install Python packages.
# This step is done before copying the rest of the code to leverage Docker's build cache.
# If requirements.txt doesn't change, this layer won't be rebuilt.
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of your application code into the container
COPY . .

# Command to run the application using Gunicorn.
# Cloud Run automatically sets the PORT environment variable to the port your container
# should listen on (typically 8080). Gunicorn is configured to bind to 0.0.0.0
# and use the value of the PORT environment variable.
# `app:app` refers to the Flask application instance named `app` within the `app.py` file.
CMD ["gunicorn", "--bind", "0.0.0.0:5000", "app:app"]
