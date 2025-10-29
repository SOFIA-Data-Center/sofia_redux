# Licensed under a 3-clause BSD style license - see LICENSE.rst

from PyQt6 import QtWidgets
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as Canvas
from matplotlib.figure import Figure
import matplotlib

matplotlib.use('QTAgg')


class MplCanvas(Canvas):
    """
    Class for handling embedding a matplotlib object in a widget.
    """
    def __init__(self):
        self.fig = Figure(figsize=(30, 35), tight_layout=True)
        Canvas.__init__(self, self.fig)
        Canvas.setSizePolicy(self,
                             QtWidgets.QSizePolicy.Policy.Expanding,
                             QtWidgets.QSizePolicy.Policy.Expanding)
        Canvas.updateGeometry(self)


class MplWidget(QtWidgets.QFrame):
    """
    Class for handling widget wrapper for plots
    """
    def __init__(self, parent=None):
        QtWidgets.QFrame.__init__(self, parent)
        self.setFrameShape(QtWidgets.QFrame.Shape.StyledPanel)

        self.canvas = MplCanvas()
        self.vbl = QtWidgets.QVBoxLayout()
        self.vbl.setContentsMargins(0, 0, 0, 0)
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)
