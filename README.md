# GenBank to GTF/GFF3 Converter

A robust web-based tool for converting GenBank flat files (`.gb`, `.gbk`) into standard Gene Transfer Format (GTF) or General Feature Format Version 3 (GFF3). This application provides a convenient API for programmatic conversion and is containerized for easy deployment.

## 🧬 Bioinformatics Context

GenBank is a comprehensive public database of nucleotide sequences and their protein translations. While rich in information, its flat-file format is often not directly compatible with many bioinformatics tools used for downstream analysis, such as genome browsers, RNA-seq quantification pipelines, or variant annotation tools. These tools frequently require annotation files in standardized formats like GTF or GFF3.

* **GTF (Gene Transfer Format)**: Primarily used for gene annotation, especially in RNA-seq analysis. It defines features like genes, transcripts, exons, and CDS, and is often used by tools like Cufflinks, StringTie, and featureCounts.

* **GFF3 (General Feature Format Version 3)**: A more flexible and hierarchical format for describing genomic features. It supports parent-child relationships between features (e.g., a gene containing multiple mRNAs, which in turn contain exons and CDS). GFF3 is widely used by genome browsers (e.g., IGV, JBrowse) and for general genomic annotation.

This converter bridges the gap between GenBank's detailed but less structured format and the structured requirements of GTF/GFF3, enabling seamless integration of GenBank annotations into various bioinformatics workflows. It meticulously handles feature locations, strands, and extracts relevant qualifiers to populate the attributes required by GTF and GFF3 standards.

## ✨ Features

* **GenBank to GTF Conversion**: Converts GenBank features (`gene`, `mRNA`, `CDS`, `exon`, `tRNA`, `rRNA`) into GTF format, including `gene_id`, `transcript_id`, `gene_name`, `product`, and `protein_id` attributes.

* **GenBank to GFF3 Conversion**: Converts GenBank features into GFF3 format, establishing hierarchical relationships (e.g., `gene` as parent of `mRNA`, `mRNA` as parent of `exon` and `CDS`). Includes `ID`, `Parent`, `Name`, `Dbxref`, `Note`, `gbkey`, `gene_biotype`, `pseudo`, `exception`, and `transl_table` attributes.

* **Robust Qualifier Extraction**: Extracts a comprehensive set of qualifiers from GenBank features to enrich the output GTF/GFF3 attributes.

* **Accurate Partiality Handling**: Correctly identifies and annotates 5' and 3' partial features in GFF3 using Biopython's precise position types.

* **Dynamic ID Generation**: Generates consistent and unique IDs for features, including fallback IDs when standard qualifiers are absent.

* **Web API**: Provides a RESTful API endpoint for file uploads or raw data paste, making it suitable for integration into other web applications or scripts.

* **Temporary File Management**: Automatically cleans up converted files from temporary storage after a configurable lifetime.

* **Containerized Deployment**: Ready for deployment using Docker, making it portable and scalable on platforms like Google Cloud Run.

## 🚀 Getting Started (Local Development)

Follow these steps to set up and run the converter locally for development or testing.

### Prerequisites

* Python 3.9+

* `pip` (Python package installer)

* `git` (for cloning the repository)

* `docker` (optional, for containerized local testing)

### Installation

1.  **Clone the repository:**

    ```bash
    git clone [https://github.com/mahdimorshediyekta/GenBank-to-GTF-GFF3-Converter.git](https://github.com/mahdimorshediyekta/GenBank-to-GTF-GFF3-Converter.git)
    cd GenBank-to-GTF-GFF3-Converter
    ```

    (Note: Replace `mahdimorshediyekta/GenBank-to-GTF-GFF3-Converter.git` with your actual repository URL if it's different).

2.  **Create a virtual environment (recommended):**

    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows: venv\Scripts\activate
    ```

3.  **Install dependencies:**

    ```bash
    pip install -r requirements.txt
    ```

### Running the Flask Application Locally

After installation, you can run the Flask development server:

```bash
python test.py
The application will start on http://0.0.0.0:5000. You can test the API by making requests to http://localhost:5000/api/convert.

💻 API Usage
The converter exposes a single primary API endpoint for conversion: /api/convert.

POST /api/convert
This endpoint accepts GenBank data either as a file upload or raw text and converts it to the specified format.

Request Body:

You can send data in one of two ways:

Multipart Form Data (for file upload):

file: The GenBank file (.gb, .gbk, .genbank).

output_format: (Optional) The desired output format. Can be gtf or gff3. Defaults to gff3.

Form Data (for raw text paste):

genbank_data: A string containing the raw GenBank file content.

output_format: (Optional) The desired output format. Can be gtf or gff3. Defaults to gff3.

Response:

A JSON object containing:

success: true if conversion was successful, false otherwise.

converted_data: The converted annotation data as a string (for smaller files).

download_url: A URL from which the converted file can be downloaded.

error: An error message if success is false.

Example curl Command (File Upload):

Bash

curl -X POST \
  -F "file=@/path/to/your/input.gb" \
  -F "output_format=gtf" \
  http://localhost:5000/api/convert
Example curl Command (Raw Data):

Bash

curl -X POST \
  -d "genbank_data=LOCUS       NC_000913             4639675 bp    DNA     circular BCT 01-JUL-2008\nDEFINITION  Escherichia coli K-12 MG1655 complete genome.\n..." \
  -d "output_format=gff3" \
  http://localhost:5000/api/convert
(Replace ... with actual GenBank content)

GET /api/download/<filename>
This endpoint serves the converted file for download. The download_url provided by the /api/convert endpoint will point to this.

🐳 Deployment
This application is designed for easy containerization and deployment.

Building the Docker Image
Ensure you have Docker installed. Navigate to the root directory of the project (where Dockerfile is located) and run:

Bash

docker build -t genbank-converter .
Running Docker Locally
To run the built Docker image on your local machine:

Bash

docker run -p 5000:5000 genbank-converter
This will make the application accessible at http://localhost:5000.

Deploying to Google Cloud Run
The provided Dockerfile is optimized for deployment on Google Cloud Run. Cloud Run automatically sets the PORT environment variable (typically to 8080), and the Gunicorn CMD in the Dockerfile is configured to listen on this variable.

Ensure gcloud CLI is configured for your Google Cloud project.

Deploy using gcloud run deploy:

Bash

gcloud run deploy genbank-to-gtf-gff3-converter \
  --source . \
  --region YOUR_REGION \
  --allow-unauthenticated # Or configure authentication as needed
Replace YOUR_REGION with your desired Google Cloud region (e.g., us-central1, europe-west1).

Cloud Run will automatically build the Docker image from your source code and deploy it.

📁 Project Structure
test.py: The main Flask application containing the GenBank conversion logic and API endpoints.

requirements.txt: Lists all Python dependencies required by the application.

Dockerfile: Instructions for building the Docker image for containerization.

temp_files/: A directory created by the application to temporarily store converted files before download. This directory is automatically managed by a background cleanup scheduler.

🤝 Contributing
Contributions are welcome! If you find a bug or have a feature request, please open an issue. If you'd like to contribute code, please fork the repository and submit a pull request.

📄 License
This project is licensed under the MIT License - see the LICENSE file for details.
