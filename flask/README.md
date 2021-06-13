# VXRD: Virtual X-Ray Diffraction
VXRD is a web interface to visualize PyXtal's [PXRD](https://nbviewer.jupyter.org/github/qzhu2017/PyXtal/blob/master/examples/tutorials_notebook/03_pxrd.ipynb) module for X-ray diffraction analysis.  The web app is currently hosted at https://vxrd.physics.unlv.edu/ for general use.

To run and view the VXRD web app on your local machine, ensure PyXtal is installed (`pip install pyxtal`) and run the following shell commands:
```bash
$ cd PyXtal/flask/
$ pip install -r requirements.txt
$ flask run
```

Add the following to your .bashrc (or equivalent)
```
# VXRD Flask Environment Variables
export SECRET_KEY=asdf1234
```

If everything is setup correctly, you should see the following output:
```bash
 * Serving Flask app "vxrd.py"
 * ...
 * Running on http://127.0.0.1:5000/ (Press CTRL+C to quit)
```

Then, open your web browser and enter the following URL:
`http://localhost:5000/`

When finished, press `CTRL+C` in your terminal to shutdown the web-app.

## VXRD: JSmol
In order to see the 3D structure visualized with JSmol, you'll need to unzip `jsmol.zip` into the following directory:
```bash
$ unzip jsmol.zip ./flask/app/static/
```
