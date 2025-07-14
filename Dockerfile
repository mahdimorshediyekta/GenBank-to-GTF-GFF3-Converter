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

# Copy the requirements file and install Python packages 
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt 

# Copy the rest of your application code into the container 
COPY . .

# Command to run the application using Gunicorn 
# This now uses the PORT environment variable, which is standard for Cloud Run and other platforms.
CMD exec gunicorn --bind 0.0.0.0:$PORT  --timeout 0 app:app
