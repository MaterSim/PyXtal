import io
import os
import plotly.graph_objects as go
from flask import render_template, flash, session, Markup
from app import app
from app.forms import MainForm
from werkzeug.utils import secure_filename
from pyxtal.XRD import XRD
from pyxtal.XRD import Similarity
from pyxtal import pyxtal
#from ase.io import read

@app.route('/')
@app.route('/index')
def index():
    return render_template('index.html')

@app.route('/individual', methods=['GET', 'POST'])
def individual():
    form = MainForm()
    if form.validate_on_submit():
        if form.upload.data:
            process_upload(form)
        if not session.get("SAVEPATH"): # new session
            flash(Markup('<span class="glyphicon glyphicon-exclamation-sign"\
                aria-hidden="true"></span><span class="sr-only">Error:</span>\
                Please upload an input file.'), 'danger alert-dismissible')
            return render_template('individual.html',
                title='Individual',
                form=form)
        else:
            process_form(form)
            return render_template('individual.html', 
                title='Individual',
                form=form,
                jsmol=True,
                plot=plot())
    # initial page visit
    return render_template('individual.html',
        title='Individual',
        form=form)

@app.route('/comparison', methods=['GET', 'POST'])
def comparison():
    form = MainForm()
    if form.validate_on_submit():
        if form.upload.data:
            process_upload(form)
        if form.upload2.data:
            process_upload(form, True)
        if not session.get("SAVEPATH")\
            or not session.get("SAVEPATH2"): # new session
            flash(Markup('<span class="glyphicon glyphicon-exclamation-sign"\
                aria-hidden="true"></span><span class="sr-only">Error:</span>\
                Please upload <strong>two (2)</strong> input files.'),
                'danger alert-dismissible')
            return render_template('comparison.html',
                title='Comparison',
                form=form)
        else:
            process_form(form, True)
            return render_template('comparison.html', 
                title='Comparison',
                form=form,
                jsmol=True,
                plot=compare())
    # initial page visit
    return render_template('comparison.html',
        title='Comparison',
        form=form)

@app.route('/struct/<type>')
def structure(type: str):
    """Return atomic structure as cif"""
    #s = pyxtal()
    #s.from_seed(session.get("SAVEPATH"))
    #struct = s.to_ase()
    struct = read(session.get("SAVEPATH"))

    if type == 'cif':
        fd = io.BytesIO()
        struct.write(fd, 'cif')
        return fd.getvalue(), 200, []
    else:
        return render_template('404.html'), 404

@app.route('/struct2/<type>')
def structure2(type: str):
    """Return 2nd atomic structure as cif"""
    struct2 = read(session.get("SAVEPATH2"))
    #s = pyxtal()
    #s.from_seed(session.get("SAVEPATH2"))
    #struct2 = s.to_ase()

    if type == 'cif':
        fd = io.BytesIO()
        struct2.write(fd, 'cif')
        return fd.getvalue(), 200, []
    else:
        1 / 0 # force error for invalid URLs

def process_upload(form, comp=False):
    """
    Save upload, check validity, and update session.
    """
    # Save uploaded file to disk
    if comp:
        f = form.upload2.data
    else:
        f = form.upload.data
    savepath = os.path.join(app.instance_path, 
        'uploads', secure_filename(f.filename))
    f.save(savepath)

    # Check readability
    try:
        #read(savepath) # attempt ase.io.read

        # Update session keys
        if comp:
            session["FILENAME2"] = f.filename
            session["SAVEPATH2"] = savepath
        else:
            session["FILENAME"] = f.filename
            session["SAVEPATH"] = savepath
            flash(Markup('<span class="glyphicon glyphicon-ok-sign"\
                aria-hidden="true"></span><span class="sr-only">Success:\
                </span> <strong>{}</strong> successfully processed.')\
                .format(session.get("FILENAME")), 
                'success alert-dismissible')
    except:
        flash(Markup('<span class="glyphicon glyphicon-exclamation-sign"\
            aria-hidden="true"></span><span class="sr-only">Error:</span> \
            Unable to read <strong>{}</strong>. Please try again or a \
            different file.').format(f.filename), 'danger alert-dismissible')

def process_form(form, comp=True):
    """
    Advanced form validation and session update.
    """
    # Retrieve form data
    max2theta = form.max2theta.data
    min2theta = form.min2theta.data

    # Advanced form-level validation
    if min2theta > max2theta:
        min2theta = 0 # use default
        flash(Markup('<span class="glyphicon glyphicon-warning-sign"\
            aria-hidden="true"></span><span class="sr-only">Error:</span> \
            2<i>&theta;</i><sub>min</sub> <strong>greater</strong> than\
            2<i>&theta;</i><sub>max</sub>&mdash;defaulting\
            2<i>&theta;</i><sub>min</sub> to 0&deg;.'), 'warning alert-dismissible')

    # Update session keys
    session["WAVELENGTH"] = form.wavelength.data
    session["MIN2THETA"] = min2theta
    session["MAX2THETA"] = max2theta
    session["RES"] = form.res.data
    session["METHOD"] = form.method.data
    session["FWHM"] = form.fwhm.data
    session["U"] = form.u.data
    session["V"] = form.v.data
    session["W"] = form.w.data
    session["A"] = form.a.data
    session["ETA_H"] = form.eta_h.data
    session["ETA_L"] = form.eta_l.data
    if comp:
        session["SHIFT"] = form.shift.data
        
def plot():
    """
    Process and return PXRD plotly.
    """
    method = session.get("METHOD")
    if method == 'gaussian' or method == 'lorentzian'\
        or method == 'pseudo-voigt':
        kwargs = {
                    'FWHM': session.get("FWHM")
                }
    elif method == 'mod_pseudo-voigt':
        kwargs = {
                    'U': session.get("U"), 
                    'V': session.get("V"),
                    'W': session.get("W"),
                    'A': session.get("A"),
                    'eta_h': session.get("ETA_H"),
                    'eta_l': session.get("ETA_L"),
                }
    struct = read(session.get("SAVEPATH"))
    xrd = XRD(struct,
        wavelength=session.get("WAVELENGTH"),
        thetas=[session.get("MIN2THETA"),
            session.get("MAX2THETA")]) 
    xrd.get_profile(method=method,
        res=session.get("RES"),
        user_kwargs=kwargs)
    flash(Markup('<span class="glyphicon glyphicon-info-sign"\
            aria-hidden="true"></span><span class="sr-only">Error:</span> \
            Showing <strong>{}</strong> with <i>{}</i> profiling.').format(
            session.get("FILENAME"), method), 'info')
    return xrd.plotly_pxrd(profile=method, height=450)

def compare():
    """
    Process and return comparison PXRD plotly.
    """
    method = session.get("METHOD")
    if method == 'gaussian' or method == 'lorentzian'\
        or method == 'pseudo-voigt':
        kwargs = {
                    'FWHM': session.get("FWHM")
                }
    elif method == 'mod_pseudo-voigt':
        kwargs = {
                    'U': session.get("U"), 
                    'V': session.get("V"),
                    'W': session.get("W"),
                    'A': session.get("A"),
                    'eta_h': session.get("ETA_H"),
                    'eta_l': session.get("ETA_L"),
                }
    files = [session.get("FILENAME"),
            session.get("FILENAME2")]
    structs = [read(session.get("SAVEPATH")),
                read(session.get("SAVEPATH2"))]
    xrds = []

    for struct in structs:
        xrd = XRD(struct,
                wavelength=session.get("WAVELENGTH"),
                thetas=[session.get("MIN2THETA"),
                session.get("MAX2THETA")])
        xrd.get_profile(method=method,
            res=session.get("RES"),
            user_kwargs=kwargs)
        xrds.append(xrd)

    S = Similarity(xrds[0].get_profile(),
        xrds[1].get_profile(),
        l=session.get("SHIFT"))

    #S.calculate()
    title = 'PXRD Similarity {:6.3f} with shift\
        {:6.3f}'.format(S.S, S.l)
    traces = []

    #for i, xrd in enumerate(xrds):
    #    traces.append(go.Scatter(x=xrd.spectra[0], y=xrd.spectra[1], name=str(files[i])))
    traces.append(go.Scatter(x=S.fx, y=S.fy, name=str(files[0])))
    traces.append(go.Scatter(x=S.fx, y=-S.gy, name=str(files[1])))
    
    fig = go.Figure(data=traces)
    fig.update_layout(
        height=450,
        title_text = title,
        title_x=0.5,
        xaxis_title = '2&#952; ({:.4f}\
        &#8491;)'.format(session.get("WAVELENGTH")),
        yaxis_title = 'Intensity')
    flash(Markup('<span class="glyphicon glyphicon-info-sign"\
            aria-hidden="true"></span><span class="sr-only">Error:</span> \
            Comparing <strong>{}</strong> and <strong>{}</strong> with\
            <i>{}</i> profiling.').format(session.get("FILENAME"),
            session.get("FILENAME2"), method), 'info')
    return fig.to_html()

def read(savepath):
    s = pyxtal()
    s.from_seed(savepath)
    return s.to_ase()
