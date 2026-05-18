'''Common utilities for various GUI classes'''

import os
try:
    from PyQt5 import QtGui, QtWidgets # pylint: disable=import-error
    from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT #pylint: disable=no-name-in-module
except ImportError:
    import sip
    sip.setapi('QString', 2)
    from PyQt4 import QtGui
    from PyQt4 import QtGui as QtWidgets # pylint: disable=import-error
    from matplotlib.backends.backend_qt4agg import FigureCanvas, NavigationToolbar2QT #pylint: disable=no-name-in-module

def add_scroll_hbox(gui, hbox):
    '''Adds buttons to hbox for scrolling through frames
    
    Connects to the class's update functions
    '''
    button = QtWidgets.QPushButton('Prev', gui)
    button.clicked.connect(gui._prev_frame) # pylint: disable=protected-access
    hbox.addWidget(button)
    button = QtWidgets.QPushButton('Next', gui)
    button.clicked.connect(gui._next_frame) # pylint: disable=protected-access
    hbox.addWidget(button)
    button = QtWidgets.QPushButton('Random', gui)
    button.clicked.connect(gui._rand_frame) # pylint: disable=protected-access
    hbox.addWidget(button)

def add_class_hbox(gui, vbox):
    hbox = QtWidgets.QHBoxLayout()
    vbox.addLayout(hbox)
    gui.class_fname = QtWidgets.QLineEdit(gui.classes.fname, gui)
    gui.class_fname.editingFinished.connect(gui._update_name)
    hbox.addWidget(gui.class_fname)
    button = QtWidgets.QPushButton('Save Classes', gui)
    button.clicked.connect(gui.classes.save)
    hbox.addWidget(button)
    hbox.addStretch(1)

class MyNavigationToolbar(NavigationToolbar2QT):
    def _icon(self, name, color=None):
        fname = os.path.abspath(os.path.dirname(__file__) + '/icons/' + name)
        pm = QtGui.QPixmap(fname)
        #if hasattr(pm, 'setDevicePixelRatio'):
        #    pm.setDevicePixelRatio(self.canvas._dpi_ratio)
        return QtGui.QIcon(pm)
