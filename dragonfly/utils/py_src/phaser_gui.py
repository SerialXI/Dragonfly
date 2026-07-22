import os
try:
    from PyQt5 import QtCore, QtWidgets, QtGui # pylint: disable=import-error
    from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT #pylint: disable=no-name-in-module
    os.environ['QT_API'] = 'pyqt5'
except ImportError:
    import sip
    sip.setapi('QString', 2)
    from PyQt4 import QtCore, QtGui # pylint: disable=import-error
    from PyQt4 import QtGui as QtWidgets # pylint: disable=import-error
    from matplotlib.backends.backend_qt4agg import FigureCanvas, NavigationToolbar2QT #pylint: disable=no-name-in-module
    os.environ['QT_API'] = 'pyqt'

import os.path as op
import numpy as np
import h5py
from matplotlib import colors
from matplotlib.figure import Figure

from . import gui_utils
from . import class_phaser


class PhaserWorker(QtCore.QObject):
    progress = QtCore.pyqtSignal(str)
    completed = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str)
    done = QtCore.pyqtSignal()

    def __init__(self, intens, num_supp, positivity, algorithms, num_runs):
        super().__init__()
        self.intens = intens
        self.num_supp = num_supp
        self.positivity = positivity
        self.algorithms = algorithms
        self.num_runs = num_runs

    @QtCore.pyqtSlot()
    def run(self):
        try:
            phaser = class_phaser.ClassPhaser(self.intens, num_supp=self.num_supp,
                                              positivity=self.positivity)
            phaser.phase(self.algorithms, num_runs=self.num_runs,
                         progress_callback=self.progress.emit)
        except Exception as exc: # pylint: disable=broad-except
            self.failed.emit('%s: %s' % (type(exc).__name__, exc))
        else:
            self.completed.emit(phaser)
        finally:
            self.done.emit()


class Phaser2D(QtWidgets.QMainWindow):
    windowClosed = QtCore.pyqtSignal()

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.output_fname = parent.fname.text()
        self.intens = parent.vol_plotter.vol
        self.curr_intens = None
        self.phaser = None
        self._phase_thread = None
        self._phase_worker = None
        self.preprocessed = False
        self.parent.vol_plotter._get_intrad()
        self.intrad = self.parent.vol_plotter.intrad

        self._init_ui()

    def _init_ui(self):
        if self.parent.css is not None:
            self.setStyleSheet(self.parent.css)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle('Mode Phaser')
        self.window = QtWidgets.QWidget()
        vbox = QtWidgets.QVBoxLayout()
        self.window.setLayout(vbox)
        self.setCentralWidget(self.window)
        self.window.setObjectName('frame')

        self.fig = Figure(figsize=(6, 6))
        self.fig.subplots_adjust(left=0.05, right=0.99, top=0.9, bottom=0.05)
        self.canvas = FigureCanvas(self.fig)
        self.navbar = gui_utils.MyNavigationToolbar(self.canvas, self)
        vbox.addWidget(self.navbar)
        vbox.addWidget(self.canvas, stretch=1)

        line = QtWidgets.QHBoxLayout()
        vbox.addLayout(line)
        out_fname = op.basename(op.dirname(self.output_fname)) + '/' + op.basename(self.output_fname)
        label = QtWidgets.QLabel('Output file: %s'%out_fname, self)
        line.addWidget(label)
        label = QtWidgets.QLabel('(%d 2D averages)'%self.intens.shape[0], self)
        line.addWidget(label)
        line.addStretch(1)
        label = QtWidgets.QLabel('Class ', self)
        line.addWidget(label)
        self.class_num = QtWidgets.QSpinBox(self)
        self.class_num.setMinimum(0)
        self.class_num.setMaximum(self.intens.shape[0]-1)
        self.class_num.setValue(self.parent.modenum.value())
        self.class_num.valueChanged.connect(self._class_num_changed)
        line.addWidget(self.class_num)

        line = QtWidgets.QHBoxLayout()
        vbox.addLayout(line)
        label = QtWidgets.QLabel('Radius Min:', self)
        line.addWidget(label)
        self.radmin = QtWidgets.QLineEdit(
            str(self.parent.settings.value('phaser2d/radmin', defaultValue='15')), self)
        self.radmin.setFixedWidth(30)
        self.radmin.returnPressed.connect(self._preprocess)
        line.addWidget(self.radmin)
        label = QtWidgets.QLabel('Max:', self)
        line.addWidget(label)
        self.radmax = QtWidgets.QLineEdit(
            str(self.parent.settings.value('phaser2d/radmax',
                                           defaultValue=str(self.intens.shape[-1]//2-1))), self)
        self.radmax.setFixedWidth(30)
        self.radmax.returnPressed.connect(self._preprocess)
        line.addWidget(self.radmax)
        label = QtWidgets.QLabel('Kernel:', self)
        line.addWidget(label)
        self.kwidth = QtWidgets.QLineEdit(
            str(self.parent.settings.value('phaser2d/kwidth', defaultValue='15')), self)
        self.kwidth.setFixedWidth(30)
        self.kwidth.returnPressed.connect(self._preprocess)
        line.addWidget(self.kwidth)
        line.addStretch(1)
        self.preprocess_button = QtWidgets.QPushButton('Preprocess', self)
        self.preprocess_button.clicked.connect(self._preprocess)
        line.addWidget(self.preprocess_button)

        line = QtWidgets.QHBoxLayout()
        vbox.addLayout(line)
        label = QtWidgets.QLabel('Num support:', self)
        line.addWidget(label)
        self.num_supp = QtWidgets.QLineEdit(
            str(self.parent.settings.value('phaser2d/num_supp', defaultValue='1000')), self)
        self.num_supp.setFixedWidth(60)
        line.addWidget(self.num_supp)
        label = QtWidgets.QLabel('Algorithm:', self)
        line.addWidget(label)
        self.algo_str = QtWidgets.QLineEdit(
            str(self.parent.settings.value('phaser2d/algo_str', defaultValue='50 ER 100 DM 100 ER')), self)
        self.algo_str.setFixedWidth(180)
        line.addWidget(self.algo_str)
        line.addStretch(1)
        self.pos_flag = QtWidgets.QCheckBox('Positivity', self)
        self.pos_flag.setChecked(self._get_bool_setting('phaser2d/positivity', True))
        line.addWidget(self.pos_flag)
        self.phase_button = QtWidgets.QPushButton('Phase', self)
        self.phase_button.clicked.connect(self._phase)
        line.addWidget(self.phase_button)

        line = QtWidgets.QHBoxLayout()
        vbox.addLayout(line)
        self.phasing_status = QtWidgets.QLabel('', self)
        line.addWidget(self.phasing_status, stretch=1)
        label = QtWidgets.QLabel('Runs:', self)
        line.addWidget(label)
        self.num_runs = QtWidgets.QLineEdit(
            str(self.parent.settings.value('phaser2d/num_runs', defaultValue='1')), self)
        self.num_runs.setFixedWidth(40)
        line.addWidget(self.num_runs)
        self.show_icalc = QtWidgets.QCheckBox('Show I_calc', self)
        self.show_icalc.stateChanged.connect(self._plot)
        self.show_icalc.setEnabled(False)
        line.addWidget(self.show_icalc)
        self.show_supp = QtWidgets.QCheckBox('Show support', self)
        self.show_supp.setChecked(self._get_bool_setting('phaser2d/show_support', False))
        self.show_supp.stateChanged.connect(self._plot)
        line.addWidget(self.show_supp)
        self.save_button = QtWidgets.QPushButton('Save', self)
        self.save_button.clicked.connect(self._save)
        line.addWidget(self.save_button)

        self._class_num_changed(self.class_num.value())
        self.show()

    def _class_num_changed(self, num):
        self.curr_intens = self.intens[num]
        self.preprocessed = False
        self.phaser = None
        self.show_icalc.setEnabled(False)
        self._plot()

    def _plot(self, state=None):
        view_limits = None
        if state is not None:
            view_limits = [(ax.get_xlim(), ax.get_ylim()) for ax in self.fig.axes]

        exponent = self.parent.expstr.text()
        rangemin = float(self.parent.rangemin.text())
        rangemax = float(self.parent.rangestr.text())
        if exponent == 'log':
            norm = colors.SymLogNorm(linthresh=rangemax*1.e-2, vmin=rangemin, vmax=rangemax)
        else:
            norm = colors.PowerNorm(float(exponent), vmin=rangemin, vmax=rangemax)
        cmap = self.parent.color_map.checkedAction().text()
        size = self.curr_intens.shape[-1]
        cen = size // 2
        if self.show_icalc.isChecked() and self.phaser is not None:
            dens = self.phaser.proj_direct(self.phaser.current)
            plot_intens = np.abs(np.fft.fftshift(np.fft.fftn(dens)))**2
        else:
            plot_intens = self.curr_intens.copy()
            plot_intens[plot_intens<0] = np.nan

        self.fig.clear()
        ax = self.fig.add_subplot(111)
        ax.imshow(plot_intens, extent=[-cen-0.5, cen+0.5, cen+0.5, -cen-0.5], norm=norm, cmap=cmap)
        ax.set_facecolor('dimgray')

        if self.phaser is not None:
            s, e = size//3, 2*size//3
            dens = self.phaser.current[s:e, s:e]
            supp = self.phaser.support[s:e, s:e]

            ax = self.fig.add_axes([0.7, 0.7, 0.29, 0.29])
            if self.show_supp.isChecked():
                alpha = np.clip(supp.astype('f8'), 0, 1)
                alpha = 0.2 + 0.8 * alpha
                alpha[supp <= 0] = 0.2
                ax.imshow(np.ones(alpha.shape), vmax=2, vmin=0)
                ax.imshow(dens, alpha=alpha, cmap='gray_r', interpolation='gaussian')
            else:
                ax.imshow(dens, cmap='gray_r', interpolation='gaussian')
            ax.set_xticks([])
            ax.set_yticks([])

        if view_limits is not None:
            for ax, (xlim, ylim) in zip(self.fig.axes, view_limits):
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)

        self.canvas.draw()

    def _preprocess(self):
        raw_intens = self.intens[self.class_num.value()]
        radmin = float(self.radmin.text())
        radmax = float(self.radmax.text())
        kwidth = float(self.kwidth.text())
        self.curr_intens = class_phaser.preproc(raw_intens, self.intrad, radmin, radmax, kwidth)
        self.preprocessed = True
        self._plot()

    def _phase(self):
        if not self.preprocessed:
            self.phasing_status.setText('Preprocess intensity first')
            return
        if self._phase_thread is not None:
            return
        self.phasing_status.setText('Starting phasing...')

        algo = self._get_algo_list()
        num_runs = int(self.num_runs.text())
        self._set_phasing_controls_enabled(False)

        self._phase_thread = QtCore.QThread(self)
        self._phase_worker = PhaserWorker(self.curr_intens.copy(), int(self.num_supp.text()),
                                          self.pos_flag.isChecked(), algo, num_runs)
        self._phase_worker.moveToThread(self._phase_thread)
        self._phase_thread.started.connect(self._phase_worker.run)
        self._phase_worker.progress.connect(self.phasing_status.setText)
        self._phase_worker.completed.connect(self._phase_completed)
        self._phase_worker.failed.connect(self._phase_failed)
        self._phase_worker.done.connect(self._phase_thread.quit)
        self._phase_worker.done.connect(self._phase_worker.deleteLater)
        self._phase_thread.finished.connect(self._phase_thread_finished)
        self._phase_thread.finished.connect(self._phase_thread.deleteLater)
        self._phase_thread.start()

    def _phase_completed(self, phaser):
        self.phaser = phaser
        self.phasing_status.setText('Phasing complete')
        self.show_icalc.setEnabled(True)
        self._plot()

    def _phase_failed(self, message):
        self.phasing_status.setText('Phasing failed: %s' % message)

    def _phase_thread_finished(self):
        self._phase_worker = None
        self._phase_thread = None
        self._set_phasing_controls_enabled(True)

    def _set_phasing_controls_enabled(self, enabled):
        controls = (self.class_num, self.radmin, self.radmax, self.kwidth,
                    self.preprocess_button, self.num_supp, self.algo_str,
                    self.pos_flag, self.phase_button, self.num_runs, self.save_button)
        for control in controls:
            control.setEnabled(enabled)

    def _get_algo_list(self):
        algo = []
        tokens = self.algo_str.text().split()
        tpos = 0
        while True:
            algo += int(tokens[tpos]) * [tokens[tpos+1]]
            tpos += 2
            if tpos >= len(tokens):
                break
        return algo

    def _save(self):
        if self.phaser is None:
            self.phasing_status.setText('Phase intensity first')
            return
        with h5py.File(self.output_fname, 'a') as f:
            if 'phasing' not in f:
                print('Adding phasing group to output file')
                f['phasing/preproc_intens'] = np.ones_like(self.intens) * np.nan
                f['phasing/dens'] = np.ones_like(self.intens) * np.nan
                f['phasing/support'] = np.zeros(self.intens.shape, dtype='f8')
            elif 'phasing/support' in f and f['phasing/support'].dtype == np.dtype('bool'):
                old_support = f['phasing/support'][:].astype('f8')
                del f['phasing/support']
                f['phasing/support'] = old_support

            num = self.class_num.value()
            print('Updating data for class', num)
            f['phasing/preproc_intens'][num] = self.curr_intens
            f['phasing/dens'][num] = self.phaser.current
            f['phasing/support'][num] = self.phaser.support

    def _save_settings(self):
        self.parent.settings.setValue('phaser2d/radmin', self.radmin.text())
        self.parent.settings.setValue('phaser2d/radmax', self.radmax.text())
        self.parent.settings.setValue('phaser2d/kwidth', self.kwidth.text())
        self.parent.settings.setValue('phaser2d/num_supp', self.num_supp.text())
        self.parent.settings.setValue('phaser2d/num_runs', self.num_runs.text())
        self.parent.settings.setValue('phaser2d/algo_str', self.algo_str.text())
        self.parent.settings.setValue('phaser2d/positivity', self.pos_flag.isChecked())
        self.parent.settings.setValue('phaser2d/show_support', self.show_supp.isChecked())

    def _get_bool_setting(self, name, default):
        value = self.parent.settings.value(name, defaultValue=default)
        if isinstance(value, str):
            return value.lower() in ('1', 'true', 'yes')
        return bool(value)

    def closeEvent(self, event):
        if self._phase_thread is not None and self._phase_thread.isRunning():
            self.phasing_status.setText('Wait for phasing to finish before closing')
            event.ignore()
            return
        self._save_settings()
        self.windowClosed.emit()
        event.accept()
