# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Run Redux reduction objects from a GUI interface."""

from pathlib import Path
import sys

from sofia_redux.pipeline.interface import Interface
from sofia_redux.pipeline.gui.main import ReduxMainWindow

try:
    from PyQt6 import QtWidgets, QtCore, QtGui
except ImportError:
    HAS_PYQT6 = False
    QtWidgets, QtCore, QtGui = None, None, None
else:
    HAS_PYQT6 = True

__all__ = ['Application', 'main']


class Application(Interface):
    """
    Graphical interface to Redux reduction objects.

    This class provides a Qt5 GUI that allows interactive parameter
    setting and reduction step running.  Intermediate data viewers
    are also supported.  Most functionality is inherited from the
    `Interface` class.

    Attributes
    ----------
    app: QApplication
        A top-level Qt widget.
    """

    def __init__(self, configuration=None):
        """
        Initialize the application, with an optional configuration.

        Parameters
        ----------
        configuration : `Configuration`, optional
            Configuration items to be used for all reductions
        """
        if not HAS_PYQT6:  # pragma: no cover
            raise ImportError('PyQt6 package is required for Redux GUI.')
        super().__init__(configuration)
        self.app = None

    def run(self):
        """Start up the application."""

        # Start application
        self.app = QtWidgets.QApplication(sys.argv)

        redux_icon_file = Path(__file__).parent / "gui/icons/redux_icon.png"
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(str(redux_icon_file)),
                       QtGui.QIcon.Mode.Normal,
                       QtGui.QIcon.State.Off)
        self.app.setWindowIcon(icon)
        self.app.setApplicationName('Redux')

        # Start a timer to allow the python interpreter to run occasionally
        # (without this, ctrl-c is swallowed by the event loop)
        timer = QtCore.QTimer()
        timer.start(200)
        timer.timeout.connect(lambda: None)

        # Start up the main window and event loop
        mw = ReduxMainWindow(self)
        mw.show()
        mw.raise_()
        sys.exit(self.app.exec())


def main():
    """Run the Redux GUI."""
    Application.tidy_log()
    app = Application()
    app.run()
