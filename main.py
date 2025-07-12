# main.py
from flask import Flask, request, send_file, jsonify
import os
import tempfile
from werkzeug.utils import secure_filename
from app import convert_genbank_to_annotation_format

app = Flask(__name__)

@app.route("/convert", methods=["POST"])
def convert():
    if 'file' not in request.files:
        return jsonify({"error": "No file part in the request"}), 400

    file = request.files['file']
    output_format = request.form.get("format", "gtf").lower()
    if output_format not in ["gtf", "gff3"]:
        return jsonify({"error": "Unsupported format. Use 'gtf' or 'gff3'."}), 400

    if file.filename == '':
        return jsonify({"error": "No selected file"}), 400

    filename = secure_filename(file.filename)

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, filename)
            output_path = os.path.join(tmpdir, f"converted.{output_format}")

            file.save(input_path)
            convert_genbank_to_annotation_format(input_path, output_path, output_format)

            return send_file(output_path, as_attachment=True)

    except Exception as e:
        return jsonify({"error": str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=8080)
