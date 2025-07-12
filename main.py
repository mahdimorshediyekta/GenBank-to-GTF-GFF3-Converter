# main.py
from flask import Flask, request, send_file, jsonify
from flask_cors import CORS 
import os
import tempfile
from converter_logic import convert_genbank_to_annotation_format
import logging

app = Flask(__name__)
CORS(app) # Enable CORS for all origins, or specify origins: CORS(app, resources={r"/convert": {"origins": "https://sciencecodons.com"}})
logging.basicConfig(level=logging.INFO) # Ensure logging is configured for Flask

@app.route('/convert', methods=['POST'])
def convert_file():
    if 'file' not in request.files:
        return jsonify({"error": "No file part in the request"}), 400
    file = request.files['file']
    if file.filename == '':
        return jsonify({"error": "No selected file"}), 400

    output_format = request.form.get('format', 'gtf').lower()
    if output_format not in ['gtf', 'gff3']:
        return jsonify({"error": "Invalid output format. Choose 'gtf' or 'gff3'."}), 400

    if file:
        try:
            # Use temporary files for input and output
            with tempfile.NamedTemporaryFile(delete=False, suffix=".gb") as tmp_input:
                file.save(tmp_input.name)
                input_path = tmp_input.name

            output_suffix = ".gff" if output_format == "gff3" else ".gtf"
            with tempfile.NamedTemporaryFile(delete=False, suffix=output_suffix) as tmp_output:
                output_path = tmp_output.name

            logging.info(f"Received file: {file.filename}, converting to {output_format}")
            convert_genbank_to_annotation_format(input_path, output_path, output_format)
            logging.info(f"Conversion complete. Sending {output_path}")

            return send_file(output_path, as_attachment=True, download_name=f"{os.path.splitext(file.filename)[0]}.{output_suffix.lstrip('.')}")

        except Exception as e:
            logging.error(f"Conversion error: {e}", exc_info=True)
            return jsonify({"error": f"An error occurred during conversion: {str(e)}"}), 500
        finally:
            # Clean up temporary files
            if 'input_path' in locals() and os.path.exists(input_path):
                os.unlink(input_path)
            if 'output_path' in locals() and os.path.exists(output_path):
                os.unlink(output_path)

@app.route('/')
def index():
    return "GenBank to Annotation Converter API is running!"

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=os.environ.get('PORT', 8080))
