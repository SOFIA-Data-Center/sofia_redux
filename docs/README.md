# SOFIA Redux documentation

The online documentation is available at
https://redux.sofiadatacenter.de


## Build docs

### Build environment

The documentation uses Sphinx to build HTML and PDF versions of the Redux documentation as well
as individual instrument manuals.

We use tox to manage the build environments with all requirements including Sphinx. Install tox
from your system packages, or into a dedicated virtual environment

    python3 -m venv .toxenv
    source .toxenv/bin/activate
    python -m pip install tox

See also the GitHub Actions deployment in
[`../.github/workflows/docs.yml`](../.github/workflows/docs.yml).

### HTML docs

To build all HTML docs, run:

    tox -e build-docs

### Redux User's manual PDFs

To make the manual PDFs:

    tox -e build-docs-pdfusermanual-exes
    tox -e build-docs-pdfusermanual-fifils
    tox -e build-docs-pdfusermanual-flitecam
    tox -e build-docs-pdfusermanual-forcast
    tox -e build-docs-pdfusermanual-hawc

The SDC publishes these manuals with a SDC-MAN document number and a cover sheet:
[SDC Redux User's Manuals](https://redux.sofiadatacenter.de/stable/sofia_redux/index.html#manuals)
