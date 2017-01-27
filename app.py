

#!/usr/bin/python
#export FLASK_DEBUG=1

from flask import Flask, render_template, request, redirect, url_for, send_from_directory
from werkzeug import secure_filename
import module
import subprocess
import os


app = Flask(__name__)
#app = Flask(static_folder='/Users/james/Desktop/FlaskApp/static')

# This is the path to the upload directory
app.config['UPLOAD_FOLDER'] = 'uploads/'


# These are the extension that we are accepting to be uploaded
app.config['ALLOWED_EXTENSIONS'] = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif'])

# For a given file, return whether it's an allowed type or not
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

# This route will show a form to perform an AJAX request
@app.route("/")
def main():
    return render_template('index2.html')

@app.route('/choose', methods=['POST'])
def choose():
    return render_template('choose_uploads.html')


@app.route('/upload', methods=['POST'])
def upload():
    # Get the name of the uploaded files
    uploaded_files = request.files.getlist("file[]")
    filenames = []
    g1 = request.form.get('Group 1')
    g2 = request.form.get('Group 2')
    g3 = request.form.get('Group 3')
    g4 = request.form.get('Group 4')
    groups = [g1, g2 , g3 , g4]

    for file in uploaded_files:
        # Check if the file is one of the allowed types/extensions
        if file and allowed_file(file.filename):
            # Make the filename safe, remove unsupported chars
            filename = secure_filename(file.filename)
            # Move the file form the temporal folder to the upload
            # folder we setup
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            # Save the filename into a list, we'll use it later
            filenames.append(filename)
            # Redirect the user to the uploaded_file route, which
            # will basicaly show on the browser the uploaded file
            # Load an html page with a link to each uploaded file
    return render_template('upload.html', filenames=filenames, groups=groups)



@app.route('/results', methods=['POST'])
def analyze():
    # groups = request.form.getlist('experimental group')
    # path = './uploads'
    # raw_files = os.listdir(path)[1:len(os.listdir(path))]
    # module.data_prep(raw_files)
    # blast_results = module.run_blast(raw_files)
    # module.format_blast(blast_results)
    # module.genelist()
    # module.R_analysis(groups)
    return render_template('results2.html')
    

# This route is expecting a parameter containing the name of a file. Then it will locate that 
# file on the upload directory and show it on the browser
@app.route('/uploads/<filename>')
def uploaded_file(filename=None):
    return send_from_directory(app.config['UPLOAD_FOLDER'],
                               filename)

if __name__ == "__main__":
    app.run()