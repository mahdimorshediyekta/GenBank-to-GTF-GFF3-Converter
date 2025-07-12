# Use the official lightweight Python image.
FROM python:3.9-slim-buster

# Allow statements and log messages to immediately appear in the logs
ENV PYTHONUNBUFFERED True

# Copy local code to the container image.
WORKDIR /app
COPY . .

# Install production dependencies.
RUN pip install -r requirements.txt

# Run the web service on container startup.
CMD exec gunicorn --bind 0.0.0.0:$PORT  --timeout 0 main:app
