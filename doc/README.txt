This is the documentation folder for PyXtal. We use Sphinx to automatically generate documentation.

The following command will update the rst files to include any new python files:
`sphinx-apidoc -o . ../pyxtal -e -M -f`
Once the rst files are satisfactory, or if no new Python files have been added/deleted, you can build the html by running:
`make html`
The html files will be output to /doc/_build/html .
