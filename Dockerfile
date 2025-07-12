# Use an official Python runtime as a parent image
FROM python:3.9-slim-buster

# Set the working directory in the container
WORKDIR /app

# Install Biopython
RUN pip install biopython gunicorn flask flask-cors # Ensure flask-cors is here

# Copy the current directory contents into the container at /app
COPY . /app

# Expose the port that the application will listen on
# Cloud Run expects the application to listen on the port specified by the PORT environment variable
ENV PORT 8080

# Command to run the application (using Gunicorn for a production web server)
# You'll need a small Flask/FastAPI wrapper around your existing logic
CMD exec gunicorn --bind :$PORT --timeout 0 main:app
