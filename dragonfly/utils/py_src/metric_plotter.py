import os

import h5py
import numpy as np
from PyQt5 import QtCore, QtWidgets # pylint: disable=import-error
from PyQt5 import QtGui # pylint: disable=import-error
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvas # pylint: disable=no-name-in-module

from . import gui_utils


class MetricComboBox(QtWidgets.QComboBox):
    def __init__(self, parent=None):
        super(MetricComboBox, self).__init__(parent)
        self.setView(QtWidgets.QTreeView(self))
        self.view().setHeaderHidden(True)
        self.view().setItemsExpandable(True)
        self.view().setRootIsDecorated(True)
        self._model = QtGui.QStandardItemModel(self)
        self.setModel(self._model)

    def set_metrics(self, names, include_none=False):
        self._model.clear()
        if include_none:
            self._add_item('None')
        occ_parent = None
        for name in names:
            if name.startswith('occupancy_'):
                if occ_parent is None:
                    occ_parent = QtGui.QStandardItem('Occupancies')
                    occ_parent.setSelectable(False)
                    self._model.appendRow(occ_parent)
                item = QtGui.QStandardItem(name.replace('occupancy_', 'mode '))
                item.setData(name, QtCore.Qt.UserRole)
                occ_parent.appendRow(item)
            else:
                self._add_item(name)
        self.view().expandAll()

    def currentText(self): # pylint: disable=invalid-name
        index = self.view().currentIndex()
        if index.isValid():
            item = self._model.itemFromIndex(index)
            value = item.data(QtCore.Qt.UserRole)
            if value is not None:
                return value
        return super(MetricComboBox, self).currentText()

    def findText(self, text, flags=QtCore.Qt.MatchExactly): # pylint: disable=invalid-name,unused-argument
        return self._find_index(text)

    def setCurrentIndex(self, index): # pylint: disable=invalid-name
        if isinstance(index, QtCore.QModelIndex):
            self.view().setCurrentIndex(index)
            self.setRootModelIndex(index.parent())
            super(MetricComboBox, self).setCurrentIndex(index.row())
            self.setRootModelIndex(QtCore.QModelIndex())
        else:
            super(MetricComboBox, self).setCurrentIndex(index)

    def _add_item(self, name):
        item = QtGui.QStandardItem(name)
        item.setData(name, QtCore.Qt.UserRole)
        self._model.appendRow(item)

    def _find_index(self, text):
        for row in range(self._model.rowCount()):
            item = self._model.item(row)
            value = item.data(QtCore.Qt.UserRole)
            if value == text:
                return self._model.indexFromItem(item)
            for child_row in range(item.rowCount()):
                child = item.child(child_row)
                if child.data(QtCore.Qt.UserRole) == text:
                    return self._model.indexFromItem(child)
        return -1


class MetricPlotter(QtWidgets.QMainWindow):
    windowClosed = QtCore.pyqtSignal()

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.output_fname = parent.fname.text()
        self.metrics = {}
        self.num_frames = 0
        self._settings_restored = False

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
        self.x_metric = MetricComboBox(self)
        self.x_metric.currentIndexChanged.connect(self._plot)
        line.addWidget(self.x_metric)
        line.addWidget(QtWidgets.QLabel('Y:', self))
        self.y_metric = MetricComboBox(self)
        self.y_metric.currentIndexChanged.connect(self._plot)
        line.addWidget(self.y_metric)
        line.addWidget(QtWidgets.QLabel('Color:', self))
        self.color_metric = MetricComboBox(self)
        self.color_metric.currentIndexChanged.connect(self._plot)
        line.addWidget(self.color_metric)
        line.addStretch(1)

        line = QtWidgets.QHBoxLayout()
        vbox.addLayout(line)
        line.addWidget(QtWidgets.QLabel('Plot:', self))
        self.plot_type = QtWidgets.QComboBox(self)
        self.plot_type.addItems(['Scatter', '2D histogram', 'Hexbin'])
        self.plot_type.currentIndexChanged.connect(self._plot)
        line.addWidget(self.plot_type)
        line.addWidget(QtWidgets.QLabel('Cmap:', self))
        self.cmap = QtWidgets.QComboBox(self)
        self.cmap.addItems(['coolwarm', 'magma', 'viridis', 'cividis', 'inferno', 'plasma'])
        self._set_combo_text(self.cmap,
                             self.parent.settings.value('metric_plotter/cmap', defaultValue='coolwarm'))
        self.cmap.currentIndexChanged.connect(self._plot)
        line.addWidget(self.cmap)
        line.addWidget(QtWidgets.QLabel('Bins:', self))
        self.num_bins = QtWidgets.QSpinBox(self)
        self.num_bins.setRange(5, 1000)
        self.num_bins.setValue(int(self.parent.settings.value('metric_plotter/bins', defaultValue=100)))
        self.num_bins.valueChanged.connect(self._plot)
        line.addWidget(self.num_bins)
        self.current_mode = QtWidgets.QCheckBox('Current mode only', self)
        self.current_mode.setChecked(self._get_bool_setting('metric_plotter/current_mode', False))
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
        self.x_metric.set_metrics(names)
        self.y_metric.set_metrics(names)
        self.color_metric.set_metrics(names, include_none=True)

        self._set_combo_text(self.x_metric, old_x if old_x in names else 'frame')
        default_y = old_y if old_y in names else self._default_y_metric(names)
        self._set_combo_text(self.y_metric, default_y)
        self._set_combo_text(self.color_metric, old_c if old_c in names else 'None')
        if not self._settings_restored:
            self._restore_settings(names)
            self._settings_restored = True
        for combo in (self.x_metric, self.y_metric, self.color_metric):
            combo.blockSignals(False)

        self.current_mode.setEnabled('mode' in self.metrics)
        self.status_label.setText('%d frames, %d metrics' % (self.num_frames, len(names)))

    def _default_y_metric(self, names):
        for name in ('likelihood', 'mutual_info', 'scale', 'mode', 'orientation'):
            if name in names:
                return name
        return names[0] if names else ''

    def _set_combo_text(self, combo, text):
        index = combo.findText(text)
        if isinstance(index, QtCore.QModelIndex) and index.isValid():
            combo.setCurrentIndex(index)
        elif index >= 0:
            combo.setCurrentIndex(index)

    def _restore_settings(self, names):
        xname = self.parent.settings.value('metric_plotter/x_metric', defaultValue='frame')
        yname = self.parent.settings.value('metric_plotter/y_metric', defaultValue='')
        cname = self.parent.settings.value('metric_plotter/color_metric', defaultValue='None')
        plot_type = self.parent.settings.value('metric_plotter/plot_type', defaultValue='Scatter')

        if xname in names:
            self._set_combo_text(self.x_metric, xname)
        if yname in names:
            self._set_combo_text(self.y_metric, yname)
        if cname == 'None' or cname in names:
            self._set_combo_text(self.color_metric, cname)
        self._set_combo_text(self.plot_type, plot_type)

    def _save_settings(self):
        self.parent.settings.setValue('metric_plotter/x_metric', self.x_metric.currentText())
        self.parent.settings.setValue('metric_plotter/y_metric', self.y_metric.currentText())
        self.parent.settings.setValue('metric_plotter/color_metric', self.color_metric.currentText())
        self.parent.settings.setValue('metric_plotter/plot_type', self.plot_type.currentText())
        self.parent.settings.setValue('metric_plotter/cmap', self.cmap.currentText())
        self.parent.settings.setValue('metric_plotter/bins', self.num_bins.value())
        self.parent.settings.setValue('metric_plotter/current_mode', self.current_mode.isChecked())

    def _get_bool_setting(self, name, default):
        value = self.parent.settings.value(name, defaultValue=default)
        if isinstance(value, str):
            return value.lower() in ('1', 'true', 'yes')
        return bool(value)

    def _plot(self):
        if not self.metrics or self.x_metric.count() == 0 or self.y_metric.count() == 0:
            return

        xname = self.x_metric.currentText()
        yname = self.y_metric.currentText()
        cname = self.color_metric.currentText()
        x = self.metrics[xname].astype('f8')
        y = self.metrics[yname].astype('f8')
        keep = np.isfinite(x) & np.isfinite(y)

        if 'blacklist' in self.metrics:
            keep &= self.metrics['blacklist'] == 0
        if self.parent.blacklist is not None:
            static_keep = np.zeros_like(keep)
            nstatic = min(len(static_keep), len(self.parent.blacklist))
            static_keep[:nstatic] = self.parent.blacklist[:nstatic] == 0
            keep &= static_keep
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
        cmap = self.cmap.currentText()
        if len(x) == 0:
            ax.text(0.5, 0.5, 'No frames match filters', ha='center', va='center', transform=ax.transAxes)
        elif plot_type == 'Scatter':
            if color is None:
                ax.scatter(x, y, s=8, alpha=0.7)
            else:
                artist = ax.scatter(x, y, c=color, s=8, alpha=0.7, cmap=cmap)
                self.fig.colorbar(artist, ax=ax, label=cname)
        elif plot_type == '2D histogram':
            artist = ax.hist2d(x, y, bins=self.num_bins.value(), cmap=cmap, cmin=1)
            self.fig.colorbar(artist[3], ax=ax, label='Frame count')
        else:
            artist = ax.hexbin(x, y, gridsize=self.num_bins.value(), cmap=cmap, mincnt=1)
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
        self._save_settings()
        self.windowClosed.emit()
        event.accept()
