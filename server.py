import os
from flask import Flask, flash, request, redirect, url_for
from werkzeug.utils import secure_filename
from CombLabel import run

UPLOAD_FOLDER = './uploads'

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        if 'stock' not in request.files:
            flash('No stock file part')
            return redirect(request.url)
        stock = request.files['stock']
        if 'sequence' not in request.files:
            flash('No sequence file part')
            return redirect(request.url)
        sequence = request.files['sequence']
        if stock.filename == '':
            flash('No selected stock file')
            return redirect(request.url)
        if sequence.filename == '':
            flash('No sequence selected file')
            return redirect(request.url)
        ncs = request.form['ncs']

        # TODO: Save file and get filenames and call run function
        # run(ncs, sequenceFileName, stockFileName)

        return '''
        {}
        '''.format(ncs)
    return '''
    <!doctype html>
    <title>Upload new File</title>
    <h1>Upload new File</h1>
    <form method=post enctype=multipart/form-data>
      <input type=file name=stock>
      <input type=file name=sequence>
      <select name=ncs>
        <option value=NCD2>NCD2</option>
        <option value=2H-ND2>2H-ND2</option>
        <option value=2H-ND3>2H-ND3</option>
        <option value=ALT12>ALT12</option>
        <option value=ALT16>ALT16</option>
        <option value=NC2>NC2</option>
      <select>
      <input type=submit value=Upload>
    </form>
    '''