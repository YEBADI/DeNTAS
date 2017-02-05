#!/usr/bin/python
#export FLASK_DEBUG=1

from flask import Flask, render_template, request, redirect, url_for, send_from_directory, session
from werkzeug import secure_filename
import module
import subprocess
import os
import uuid

app = Flask(__name__)
#app = Flask(static_folder='/Users/james/Desktop/FlaskApp/static')

# This is the path to the upload directory
app.config['UPLOAD_FOLDER'] = 'uploads/'


# These are the extension that we are accepting to be uploaded
app.config['ALLOWED_EXTENSIONS'] = set(['txt', 'fasta'])

# For a given file, return whether it's an allowed type or not
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

# This route will show a form to perform an AJAX request
@app.route("/")
def main():
    return render_template('index.html')

@app.route('/choose', methods=['POST'])
def choose():
    return render_template('choose_uploads.html')

@app.route('/about')
def about():
    return render_template('about.html')

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
    #now we're making a unique user ID and saving it to the session (cookie)
    session['user'] = uuid.uuid1()
    user_folder = os.path.join(app.config['UPLOAD_FOLDER'], session['user'])

    #if path doesn't exist then make it
    if not in os.path.exists(user_folder):
        os.makedirs(user_folder)
    
    for file in uploaded_files:
        # Check if the file is one of the allowed types/extensions
        if file and allowed_file(file.filename):
            # Make the filename safe, remove unsupported chars
            filename = secure_filename(file.filename)
            # Move the file form the temporal folder to the upload
            # folder we setup
            #now each user folder is unique and files are saved respectively there
            file.save(os.path.join(user_folder, filename))
            # Save the filename into a list, we'll use it later
            filenames.append(filename)
            # Redirect the user to the uploaded_file route, which
            # will basicaly show on the browser the uploaded file
            # Load an html page with a link to each uploaded file     
    return render_template('upload.html', filenames=filenames, groups=groups)


@app.route('/results', methods=['POST'])
def analyze():
    module.apocrita_upload(session['user'])
    module.apocrita_blast()
    module.apocrita_download(session['user'])
    groups = request.form.getlist('experimental group', session['user'])
    module.R_analysis(groups, session['user'])
    module.pdf_combine(session['user'])
    module.clean_up(session['user'])
    return render_template('results.html')
    
# This route is expecting a parameter containing the name of a file. Then it will locate that 
# file on the upload directory and show it on the browser
@app.route('/uploads/<filename>')
def uploaded_file(filename=None):
    return send_from_directory(app.config['UPLOAD_FOLDER'],
                               filename)

if __name__ == "__main__":
    app.secret_key = "t3tmiU+u7VXna4M98i+G3IXqjZ5ZSmrviQgbiyy5Tcw6kdwu/uZ0gw+0pLbzRsBINFjc6FjRd9QTsbByFSmS9w=="
    app.run()
