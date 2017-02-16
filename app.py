#!/usr/bin/python
#export FLASK_DEBUG=1

from flask import Flask, render_template, request, redirect, url_for, send_from_directory
from werkzeug import secure_filename
import module
import subprocess
import os
import numpy


app = Flask(__name__)

# This is the path to the upload directory
app.config['UPLOAD_FOLDER'] = 'uploads/'

# These are the extensions that we are accepting to be uploaded
app.config['ALLOWED_EXTENSIONS'] = set(['txt', 'fasta'])

# For a given file, return whether it's an allowed type or not
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

# This is the main app route and will return the index page with information
# about the project
@app.route("/")
def main():
    return render_template('index.html')

# This is app route and will return the about page with further information
@app.route('/about')
def about():
    return render_template('about.html')

# This app route returns the first analysis page in which the user inputs the names
# of their exp.groups and chooses their raw data to upload
@app.route('/choose', methods=['POST'])
def choose():
    module.clear_previous_results()       # clear previous results before a new run
    return render_template('choose_uploads.html')


# This app route returns the second analysis page in which the user selects their model
# organism, and chooses the exp. group for each sample from a drop down box
@app.route('/upload', methods=['POST'])
def upload():
    uploaded_files = request.files.getlist("file[]") # Get the names of the uploaded files
    
    g1 = request.form.get('Group 1')
    g2 = request.form.get('Group 2')
    g3 = request.form.get('Group 3')
    g4 = request.form.get('Group 4')
    groups = [g1, g2 , g3 , g4]                 #create a list of exp groups
    filenames = []                              
    for file in uploaded_files:                  
        if file and allowed_file(file.filename):       # is file extension allowed?     
            filename = secure_filename(file.filename)  # Make safe (remove unsupported chars) 
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename)) # Save the file in uploads
            filenames.append(filename)                 # add filename to list         
    # Load html page with i) links to uploaded files ii) dropdown boxes populated by groups     
    return render_template('upload.html', filenames=filenames, groups=groups)


@app.route('/results', methods=['POST'])         #see module.py for more detailed function description
def analyze():
    module.apocrita_upload()                     # copies files from local (uploads) to remote (raw) folder
    organism = request.form.get('organism')      #from dropdown box on pervious page
    module.apocrita_blast(organism)              #runs organism specific call_blast.sh over ssh
    module.apocrita_download()                   # copies files from remote (results) to local (blast) folder
    groups = request.form.getlist('experimental group')  # each sample assigned a group
    module.R_analysis(groups)
    module.pdf_combine()                         # R output collated 
    module.clean_up()                            # removes temporary files folders
    ngroups = len(numpy.unique(groups))          
    # different results pages (to accomodate different numbers of graphs) 
    if ngroups==2:
        return render_template('results2.html')
    elif ngroups==3:
        return render_template('results3.html')
    elif ngroups==4:
        return render_template('results4.html')
    
# This route is expecting a parameter containing the name of a file. Then it will locate that 
# file on the upload directory and show it on the browser
@app.route('/uploads/<filename>')
def uploaded_file(filename=None):
    return send_from_directory(app.config['UPLOAD_FOLDER'],
                               filename)

if __name__ == "__main__":
    app.run(debug=True)      #turn debugger & reloading feature on for development