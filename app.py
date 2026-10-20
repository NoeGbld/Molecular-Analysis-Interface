# -*- coding: utf-8 -*-

import os
import sys
import json
import signal
import subprocess
import io
import csv
import logging
from threading import Lock
from flask import Flask, render_template, request, Response

from core.molecule import Molecule
from core.plot_html import plot_absorption_spectrum
from core.parsers.pubchem_parser import fill_from_pubchem
from core.parsers.astroquery_parser import fill_from_radis, fill_from_astroquery
from core.parsers.chemical_parser import fill_from_chemical
from core.parsers.nist_parser import fill_from_nist
from core.parsers.mvoc_parser_html import fill_from_mvoc
from core.parsers.sander_parser import fill_from_sander
from core.tools import install_playwright_browsers

# We disable request logs for Werkzeug / Flask
log = logging.getLogger('werkzeug')
log.setLevel(logging.ERROR)

app = Flask(__name__)

# Thread-safe variables to track and cancel stale background search tasks
search_lock = Lock()
current_search_id = 0

def execute_step(name, func, mol, **kwargs):
    # Safely executes a parser function with optional keyword arguments and returns a formatted log message
    try:
        if kwargs:
            msg = func(mol, **kwargs)
        else:
            msg = func(mol)
            
        if msg is None:
            return f"[WARN] {name}: No data found."
        return msg
    except Exception as e:
        return f"[WARN] {name} error: {str(e)}"


@app.route('/progress')
def progress():
    # Server-Sent Events (SSE) stream endpoint to process parsers sequentially
    global current_search_id
    
    molecule_name = request.args.get('molecule_name', '').strip()
    temperature = float(request.args.get('temperature', 295.15))
    pressure = float(request.args.get('pressure', 1.0))/1.01325 # bar to atm
    isotope = int(request.args.get('isotope', 1))
    fmin = float(request.args.get('fmin', 30.0))
    fmax = float(request.args.get('fmax', 800.0))
    resolution = float(request.args.get('resolution', 0.1))
    
    # Safely intercept the UI checkbox from the main thread before starting the stream generator
    plot_spectre_checked = request.args.get('plot_spectre') is not None

    # Increment and assign a unique task ID to this request session
    with search_lock:
        current_search_id += 1
        my_search_id = current_search_id

    def generate():
        if not molecule_name:
            yield f"data: {json.dumps({'status': 'error', 'log': 'Molecule name missing!'})}\n\n"
            return

        try:
            # Step 1: Initialize base target molecule
            yield f"data: {json.dumps({'percent': 5, 'log': 'Creating Molecule object...'})}\n\n"
            mol = Molecule(
                molecule_name, 
                temperature=temperature,
                pressure=pressure,
                isotope=isotope,
                fmin=fmin,
                fmax=fmax,
                resolution=resolution
            )
            
            # Check for task cancellation before moving to heavy database operations
            if current_search_id != my_search_id:
                return

            # Step 2: Query PubChem primary database
            yield f"data: {json.dumps({'percent': 15, 'log': '[PubChem] Processing...'})}\n\n"
            try:
                fill_from_pubchem(mol)
                yield f"data: {json.dumps({'percent': 20, 'log': '[PubChem] OK'})}\n\n"
            except Exception as e:
                yield f"data: {json.dumps({'status': 'error', 'log': f'[PubChem] Critical fail: {e}'})}\n\n"
                return

            # Dynamic assignment of the spectrum parsing function based on UI selection
            spectral_function = fill_from_radis if plot_spectre_checked else fill_from_astroquery
            spectral_name = "RADIS" if plot_spectre_checked else "Astroquery"

            steps = [
                (spectral_name, spectral_function, 35),
                ("ChemicalBook", fill_from_chemical, 50),
                ("NIST", fill_from_nist, 65),
                ("mVOC", fill_from_mvoc, 80),
                ("Sander", fill_from_sander, 90)
            ]

            for name, func, percent in steps:
                # Intercept loop if a newer search request has been triggered by the user
                if current_search_id != my_search_id:
                    return

                yield f"data: {json.dumps({'percent': percent, 'log': f'[{name}] Processing...'})}\n\n"
                msg = execute_step(name, func, mol)
                status = "warning" if "[WARN]" in msg else "success"

                yield f"data: {json.dumps({
                    'percent': percent,
                    'log': msg,
                    'status': status
                })}\n\n"

            # Final cancellation check before rendering templates and saving state
            if current_search_id != my_search_id:
                return

            # Step 4: Finalize spectrum graph and compile results template
            yield f"data: {json.dumps({'percent': 95, 'log': 'Generating spectrum graph...'})}\n\n"
            graph_img = plot_absorption_spectrum(mol)
            
            has_graph = graph_img is not None and len(str(graph_img).strip()) > 0
            html = app.jinja_env.get_template('results.html').render(mol=mol, graph=graph_img)
            
            yield f"data: {json.dumps({
                'percent': 100, 
                'html': html, 
                'has_graph': has_graph,
                'active_params': {
                    'temperature': mol.temperature, 
                    'pressure': mol.pressure, 
                    'isotope': mol.isotope, 
                    'fmin': mol.fmin, 
                    'fmax': mol.fmax, 
                    'resolution': mol.resolution
                }
            })}\n\n"
            
        except Exception as e:
            yield f"data: {json.dumps({'status': 'error', 'log': f'Global error: {str(e)}'})}\n\n"

    return Response(generate(), mimetype='text/event-stream')


@app.route('/download-spectrum-csv')
def download_spectrum_csv():
    # Rebuilds the dataset on-demand and builds a downloadable CSV file
    molecule_name = request.args.get('molecule_name', '').strip()
    temperature = float(request.args.get('temperature', 295.15))
    pressure = float(request.args.get('pressure', 1.0))
    isotope = int(request.args.get('isotope', 1))
    fmin = float(request.args.get('fmin', 30.0))
    fmax = float(request.args.get('fmax', 800.0))
    resolution = float(request.args.get('resolution', 0.1))

    if not molecule_name:
        return "Molecule name missing", 400

    mol = Molecule(
        molecule_name, 
        temperature=temperature,
        pressure=pressure,
        isotope=isotope,
        fmin=fmin,
        fmax=fmax,
        resolution=resolution
    )
    
    try:
        fill_from_pubchem(mol)
        fill_from_radis(mol) 
    except Exception as e:
        return f"Error extracting spectrum data: {str(e)}", 500

    if not hasattr(mol.spectre, 'nu_GHz') or mol.spectre.nu_GHz is None or \
       not hasattr(mol.spectre, 'coef') or mol.spectre.coef is None or \
       len(mol.spectre.nu_GHz) == 0:
        return "No spectral data available for this molecule and range.", 404

    output = io.StringIO()
    writer = csv.writer(output, delimiter=';')
    
    writer.writerow(['Frequency (GHz)', 'Absorption coefficient'])

    for freq, coef in zip(mol.spectre.nu_GHz, mol.spectre.coef):
        writer.writerow([round(freq, 5), coef])

    output.seek(0)
    filename = f"{molecule_name}_spectrum_{fmin}_{fmax}GHz.csv".replace(" ", "_")
    
    return Response(
        output.getvalue(),
        mimetype="text/csv",
        headers={"Content-disposition": f"attachment; filename={filename}"}
    )


@app.route('/', methods=['GET', 'POST'])
def home():
    # Standard router loading the main user interface dashboard
    params = {
        'molecule_name': '', 
        'temperature': 295.15,
        'pressure': 1.0,
        'isotope': 1,
        'fmin': 30.0,
        'fmax': 800.0,
        'resolution': 0.1
    }
    return render_template('index.html', mol=None, graph=None, params=params)


@app.route('/shutdown', methods=['POST'])
def shutdown():
    # Adaptive system shutdown sequence compatible with UNIX and Windows kernels
    try:
        # Standard execution for UNIX platforms (Linux / macOS)
        os.kill(os.getpid(), signal.SIGINT)
    except ValueError:
        # Graceful fallback execution path for Windows environments
        os.kill(os.getpid(), signal.SIGTERM)
    return "Server shut down."


if __name__ == '__main__':
    # 1. Install headless browser dependencies if required
    install_playwright_browsers()
    
    # 2. Start the local server directly
    # 'host="127.0.0.1"' forces the loop locally, bypassing strict university network profiles
    app.run(host="127.0.0.1", port=5000, debug=False)