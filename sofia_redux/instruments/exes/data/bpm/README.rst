The bad pixel mask should be a single-extension FITS file containing integer
values, where 1 indicates a good pixel, 0 indicates a bad pixel, and
2 indicates a reference pixel.

Note that the default files are included in the source repository of this
package, but not in the sdists and wheels from PyPI or GitHub.
They may be downloaded separately, if desired, from the
`GitHub repository <https://github.com/SOFIA-Data-Center/sofia_redux>`__.
Otherwise, the software will attempt to automatically download and
cache the reference file as needed.
