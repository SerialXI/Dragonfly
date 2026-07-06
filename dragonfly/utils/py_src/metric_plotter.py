import os

import h5py
import numpy as np
from PyQt5 import QtCore, QtWidgets # pylint: disable=import-error
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvas # pylint: disable=no-name-in-module

from . import gui_utils


class MetricPlotter(QtWidgets.QMainWindow):
    windowClosed = QtCore.pyqtSignal()

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.output_fname = parent.fname.text()
        self.metrics = {}
        self.num_frames = 0

        self._init_ui()
        self._load_metrics()

    def _init_ui(self):
        if self.parent.css is not None:
            self.setStyleSheet(self.parent.css)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle('Frame Metrics Plotter')
        self.window = QtWidgets.QWidget()
        vbox = QtWidgets.QVBoxLayout()
        self.window.setLayout(vbox)
        self.setCentralWidget(self.window)
        self.window.setObjectName('frame')

        self.fig = Figure(figsize=(7, 6))
        self.canvas = FigureCanvas(self.fig)
        self.navbar = gui_utils.MyNavigationToolbar(self.canvas, self)
        vbox.addWidget(self.navbar)
        vbox.addWidget(self.canvas, stretch=1)

        line = QtWidgets.QHBoxLayout()
        vbox.addLayout(line)
        self.file_label = QtWidgets.QLabel('', self)
        line.addWidget(self.file_label)
        line.addStretch(1)
        self.status_label = QtWidgets.QLabel('', self)
        line.addWidget(self.status_label)

        line = QtWidgets.QHBoxLayout()
        vbox.addLayout(line)
        line.addWidget(QtWidgets.QLabel('X:', self))
        self.x_metric = QtWidgets.QComboBox(self)
        self.x_metric.currentIndexChanged.connect(self._plot)
        line.addWidget(self.x_metric)
        line.addWidget(QtWidgets.QLabel('Y:', self))
        self.y_metric = QtWidgets.QComboBox(self)
        self.y_metric.currentIndexChanged.connect(self._plot)
        line.addWidget(self.y_metric)
        line.addWidget(QtWidgets.QLabel('Color:', self))
        self.color_metric = QtWidgets.QComboBox(self)
        self.color_metric.currentIndexChanged.connect(self._plot)
        line.addWidget(self.color_metric)

        line = QtWidgets.QHBoxLayout()
        vbox.addLayout(line)
        line.addWidget(QtWidgets.QLabel('Plot:', self))
        self.plot_type = QtWidgets.QComboBox(self)
        self.plot_type.addItems(['Scatter', '2D histogram', 'Hexbin'])
        self.plot_type.currentIndexChanged.connect(self._plot)
        line.addWidget(self.plot_type)
        line.addWidget(QtWidgets.QLabel('Bins:', self))
        self.num_bins = QtWidgets.QSpinBox(self)
        self.num_bins.setRange(5, 1000)
        self.num_bins.setValue(100)
        self.num_bins.valueChanged.connect(self._plot)
        line.addWidget(self.num_bins)
        self.skip_blacklist = QtWidgets.QCheckBox('Skip blacklisted', self)
        self.skip_blacklist.setChecked(True)
        self.skip_blacklist.stateChanged.connect(self._plot)
        line.addWidget(self.skip_blacklist)
        self.current_mode = QtWidgets.QCheckBox('Current mode only', self)
        self.current_mode.stateChanged.connect(self._plot)
        line.addWidget(self.current_mode)
        button = QtWidgets.QPushButton('Refresh', self)
        button.clicked.connect(self._load_metrics)
        line.addWidget(button)
        line.addStretch(1)

        self.show()

    def update_iteration(self):
        self._load_metrics()

    def update_mode(self):
        if self.current_mode.isChecked():
            self._plot()

    def _load_metrics(self):
        self.output_fname = self.parent.fname.text()
        self.metrics = {}
        self.num_frames = 0
        self.file_label.setText('Output file: %s' % self._short_name(self.output_fname))

        if not h5py.is_hdf5(self.output_fname):
            self._set_empty('Select an HDF5 output file')
            return

        with h5py.File(self.output_fname, 'r') as fptr:
            self.num_frames = self._get_num_frames(fptr)
            if self.num_frames == 0:
                self._set_empty('No frame-wise metrics found')
                return

            self.metrics['frame'] = np.arange(self.num_frames)
            for name, dset in fptr.items():
                if isinstance(dset, h5py.Dataset) and self._is_framewise_1d(dset):
                    self.metrics[name] = dset[:]

            if 'orientations' in self.metrics:
                rots = self.metrics['orientations'].astype('i8')
                self.metrics['orientation'] = rots.copy()
                if self.parent.num_modes > 1:
                    rotind = rots // self.parent.num_modes
                    modes = rots % self.parent.num_modes
                    modes[rots < 0] = -1
                    if self.parent.num_nonrot > 0:
                        nonrot = rotind >= self.parent.num_rot
                        modes[nonrot] = rots[nonrot] - self.parent.num_modes * (self.parent.num_rot - 1)
                    self.metrics['mode'] = modes
                    self.metrics['orientation'] = rotind

            occ_dset = fptr.get('occupancies')
            if (isinstance(occ_dset, h5py.Dataset) and len(occ_dset.shape) == 2 and
                    occ_dset.shape[0] == self.num_frames):
                occ = occ_dset[:]
                for mode in range(occ.shape[1]):
                    self.metrics['occupancy_%d' % mode] = occ[:, mode]

        self._update_selectors()
        self._plot()

    def _get_num_frames(self, fptr):
        for name in ('likelihood', 'mutual_info', 'scale', 'orientations', 'blacklist'):
            if name in fptr and len(fptr[name].shape) == 1:
                return fptr[name].shape[0]
        occ_dset = fptr.get('occupancies')
        if isinstance(occ_dset, h5py.Dataset) and len(occ_dset.shape) == 2:
            return occ_dset.shape[0]
        return 0

    def _is_framewise_1d(self, dset):
        return (len(dset.shape) == 1 and dset.shape[0] == self.num_frames and
                np.issubdtype(dset.dtype, np.number))

    def _update_selectors(self):
        old_x = self.x_metric.currentText()
        old_y = self.y_metric.currentText()
        old_c = self.color_metric.currentText()
        names = sorted(self.metrics)

        for combo in (self.x_metric, self.y_metric, self.color_metric):
            combo.blockSignals(True)
            combo.clear()
        self.x_metric.addItems(names)
        self.y_metric.addItems(names)
        self.color_metric.addItem('None')
        self.color_metric.addItems(names)

        self._set_combo_text(self.x_metric, old_x if old_x in names else 'frame')
        default_y = old_y if old_y in names else self._default_y_metric(names)
        self._set_combo_text(self.y_metric, default_y)
        self._set_combo_text(self.color_metric, old_c if old_c in names else 'None')
        for combo in (self.x_metric, self.y_metric, self.color_metric):
            combo.blockSignals(False)

        self.current_mode.setEnabled('mode' in self.metrics)
        self.skip_blacklist.setEnabled('blacklist' in self.metrics)
        self.status_label.setText('%d frames, %d metrics' % (self.num_frames, len(names)))

    def _default_y_metric(self, names):
        for name in ('likelihood', 'mutual_info', 'scale', 'mode', 'orientation'):
            if name in names:
                return name
        return names[0] if names else ''

    def _set_combo_text(self, combo, text):
        index = combo.findText(text)
        if index >= 0:
            combo.setCurrentIndex(index)

    def _plot(self):
        if not self.metrics or self.x_metric.count() == 0 or self.y_metric.count() == 0:
            return

        xname = self.x_metric.currentText()
        yname = self.y_metric.currentText()
        cname = self.color_metric.currentText()
        x = self.metrics[xname].astype('f8')
        y = self.metrics[yname].astype('f8')
        keep = np.isfinite(x) & np.isfinite(y)

        if self.skip_blacklist.isChecked() and 'blacklist' in self.metrics:
            keep &= self.metrics['blacklist'] == 0
        if self.current_mode.isChecked() and 'mode' in self.metrics:
            keep &= self.metrics['mode'] == self.parent.modenum.value()
        if cname != 'None':
            color = self.metrics[cname].astype('f8')
            keep &= np.isfinite(color)
            color = color[keep]
        else:
            color = None

        x = x[keep]
        y = y[keep]
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        ax.set_facecolor('k')

        plot_type = self.plot_type.currentText()
        if len(x) == 0:
            ax.text(0.5, 0.5, 'No frames match filters', ha='center', va='center', transform=ax.transAxes)
        elif plot_type == 'Scatter':
            if color is None:
                ax.scatter(x, y, s=8, alpha=0.7)
            else:
                artist = ax.scatter(x, y, c=color, s=8, alpha=0.7, cmap='viridis')
                self.fig.colorbar(artist, ax=ax, label=cname)
        elif plot_type == '2D histogram':
            artist = ax.hist2d(x, y, bins=self.num_bins.value(), cmap='viridis', cmin=1)
            self.fig.colorbar(artist[3], ax=ax, label='Frame count')
        else:
            artist = ax.hexbin(x, y, gridsize=self.num_bins.value(), cmap='viridis', mincnt=1)
            self.fig.colorbar(artist, ax=ax, label='Frame count')

        ax.set_xlabel(xname)
        ax.set_ylabel(yname)
        ax.set_title('%s vs %s (%d frames)' % (yname, xname, len(x)))
        self.fig.tight_layout()
        self.canvas.draw()

    def _set_empty(self, message):
        self.metrics = {}
        self.x_metric.clear()
        self.y_metric.clear()
        self.color_metric.clear()
        self.status_label.setText(message)
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        ax.text(0.5, 0.5, message, ha='center', va='center', transform=ax.transAxes)
        ax.set_axis_off()
        self.canvas.draw()

    def _short_name(self, fname):
        return os.path.basename(os.path.dirname(fname)) + '/' + os.path.basename(fname)

    def closeEvent(self, event):
        self.windowClosed.emit()
        event.accept()
