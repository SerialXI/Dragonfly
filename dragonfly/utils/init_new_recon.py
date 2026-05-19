#!/usr/bin/env python

'''Initialize a new reconstruction directory.

Creates the standard Dragonfly reconstruction folder structure and either copies
the packaged default config or generates one interactively for simulation or
experimental data workflows.
'''

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import shutil
import sys

try:
    from prompt_toolkit.application import Application
    from prompt_toolkit.filters import Condition
    from prompt_toolkit.key_binding import KeyBindings
    from prompt_toolkit.formatted_text import FormattedText
    from prompt_toolkit.layout import Layout
    from prompt_toolkit.layout.containers import DynamicContainer, Float, FloatContainer, HSplit, VSplit
    from prompt_toolkit.layout.scrollable_pane import ScrollablePane
    from prompt_toolkit.layout.menus import CompletionsMenu
    from prompt_toolkit.styles import Style as PTStyle
    from prompt_toolkit.completion import PathCompleter
    from prompt_toolkit.layout.dimension import Dimension as D
    from prompt_toolkit.widgets import Box, Button, CheckboxList, Frame, Label, RadioList, TextArea
    from rich.console import Console
    HAS_INTERACTIVE_TUI = True
except ImportError:
    HAS_INTERACTIVE_TUI = False

SIM_DEFAULTS = {
    'parameters': {
        'detd': '300',
        'lambda': '6.2',
        'detsize': '150',
        'pixsize': '0.512',
        'stoprad': '10',
        'polarization': 'x',
    },
    'make_densities': {
        'pdb_code': '4BED',
        'scatt_dir': 'aux/',
        'out_density_file': 'data/densityMap.bin',
    },
    'make_intensities': {
        'in_density_file': 'make_densities:::out_density_file',
        'out_intensity_file': 'data/intensities.bin',
    },
    'make_detector': {
        'out_detector_file': 'data/det_sim.h5',
    },
    'make_data': {
        'num_data': '30000',
        'fluence': '1e11',
        'in_detector_file': 'make_detector:::out_detector_file',
        'in_intensity_file': 'make_intensities:::out_intensity_file',
        'out_photons_file': 'data/photons.emc',
    },
    'emc': {
        'in_photons_file': 'make_data:::out_photons_file',
        'in_detector_file': 'make_detector:::out_detector_file',
        'recon_type': '3d',
        'num_div': '6',
        'output_folder': 'data/',
        'log_file': 'logs/EMC.log',
        'need_scaling': '1',
        'beta_factor': '1.0',
        'beta_schedule': '2.0 10',
    },
}

EXP_DEFAULTS = {
    'emc': {
        'recon_type': '3d',
        'output_folder': 'data/',
        'log_file': 'logs/EMC.log',
        'num_div': '10',
        'need_scaling': '1',
        'beta_factor': '1.0',
        'beta_schedule': '1.41421356 10',
    },
}

RECON_DEFAULTS_2D = {
    'num_rot': '180',
    'num_modes': '8',
}

EMPTY_CHECKBOX_VALUE = '__dragonfly_empty__'

PT_APP_STYLE = PTStyle.from_dict({
    'dialog': 'bg:#000000 #d7dce2',
    'dialog frame.label': 'bg:#000000 #d7dce2 bold',
    'dialog.body': 'bg:#000000 #d7dce2',
    'dialog shadow': 'bg:#111111',
    'frame.border': '#4d535b',
    'label': '#aeb6bf',
    'input-field': 'bg:#14181d #ffffff',
    'input-field.focused': 'bg:#1f2f45 #ffffff bold',
    'footer': 'bg:#000000 #d7dce2',
    'footer-key': '#7fd7ff bold',
    'footer-mode-path': 'bg:#243447 #d7ecff bold',
    'footer-mode-edit': 'bg:#2f4f2f #e8ffe8 bold',
    'footer-info': '#d7dce2',
    'footer-hint': '#aeb6bf',
    'footer-error': '#ff8f8f bold',
    'button': 'bg:#343a40 #d7dce2',
    'button.focused': 'bg:#7a828c #ffffff bold',
    'dialog.body button': 'bg:#343a40 #d7dce2',
    'dialog.body button.focused': 'bg:#7a828c #ffffff bold',
    'dialog.body button-arrow': 'bg:#343a40 #d7dce2',
    'dialog.body button.focused button-arrow': 'bg:#7a828c #ffffff bold',
    'radio': '#aeb6bf',
    'radio-checked': '#d7dce2',
    'radio-selected': '#7fd7ff bold',
    'checkbox': '#aeb6bf',
    'checkbox-checked': '#98fb98 bold',
    'checkbox-selected': '#98fb98 bold',
    'focused': 'bg:#264f78 #ffffff bold',
}) if HAS_INTERACTIVE_TUI else None

# Shared console and filesystem helpers.

def _name_recon_dir(tag, num):
    return '%s_%04d' % (tag, num)

def _first_available_num(tag, num, prefix):
    while op.exists(op.join(prefix, _name_recon_dir(tag, num))):
        num += 1
    return num

def _create_recon_dir(tag='recon', num=1, prefix='./'):
    recon_num = _first_available_num(tag, num, prefix)
    recon_dir = op.join(prefix, _name_recon_dir(tag, recon_num))
    logging.info('New recon directory created with name: %s', recon_dir)
    os.mkdir(recon_dir)
    os.mkdir(op.join(recon_dir, 'data'))
    os.mkdir(op.join(recon_dir, 'images'))
    os.mkdir(op.join(recon_dir, 'logs'))
    link_name = _name_recon_dir(tag, recon_num)
    if not op.exists(link_name):
        os.symlink(recon_dir, link_name)
    return recon_dir

def _setup_aux_dir(recon_dir, parent_dir, copy_aux=False):
    src = op.join(parent_dir, 'aux')
    dst = op.join(recon_dir, 'aux')
    if not op.lexists(dst):
        if copy_aux:
            shutil.copytree(src, dst, symlinks=True)
        else:
            os.symlink(src, dst)

def _copy_default_config(recon_dir, parent_dir):
    src = op.join(parent_dir, 'config.ini')
    shutil.copy(src, op.join(recon_dir, 'config.ini'))

def _ordered_experimental_emc(emc):
    ordered_emc = {}
    for key in ('in_photons_file', 'in_photons_list', 'in_detector_file', 'in_detector_list'):
        if key in emc:
            ordered_emc[key] = emc[key]
    for key in ('recon_type', 'output_folder', 'log_file', 'need_scaling'):
        if key in emc:
            ordered_emc[key] = emc[key]
    for key in ('num_div', 'num_rot', 'num_modes'):
        if key in emc:
            ordered_emc[key] = emc[key]
    for key in ('beta', 'beta_factor', 'beta_schedule'):
        if key in emc:
            ordered_emc[key] = emc[key]
    for key, value in emc.items():
        if key not in ordered_emc:
            ordered_emc[key] = value
    return ordered_emc

class FullScreenWizard(object):
    '''Full-screen prompt_toolkit wizard for reconstruction setup.'''

    STEPS = {
        'directory': 'Reconstruction Directory',
        'config_style': 'Config Style',
        'workflow': 'Workflow',
        'sim_model_source': 'Simulation Model Source',
        'sim_model_value': 'Simulation Model Input',
        'sim_geometry': 'Detector Geometry',
        'sim_data': 'Simulated Data',
        'sim_recon': 'Simulation EMC Settings',
        'sim_beta': 'Simulation Beta Schedule',
        'exp_photon_mode': 'Photon Input Mode',
        'exp_photon_input': 'Photon Input Setup',
        'exp_detector': 'Experimental Detector and Output',
        'exp_recon': 'Experimental EMC Settings',
        'exp_beta': 'Experimental Beta Schedule',
    }

    # Wizard lifecycle and layout.

    def __init__(self, args, parent_dir):
        self.args = args
        self.parent_dir = parent_dir
        self.status = 'Press Ctrl-C to exit without making any changes.'
        self.status_style = 'class:footer-hint'
        self.cancelled = False
        self.edit_mode = False
        self.current_step = 'directory'
        self.current_focus = None
        self.current_container = Box(Label('Loading...'))
        self.summary_container = Box(Label('Loading...'))
        self.footer = Label(text=self._footer_mode_text, style='class:footer')
        self.footer_hint = Label(text=self._footer_status_text, style='class:footer')
        self.state = self._initial_state()
        self.kb = self._make_key_bindings()
        self.root = FloatContainer(
            content=HSplit([
                VSplit([
                    Box(ScrollablePane(DynamicContainer(lambda: self.summary_container)), width=D(weight=1)),
                    Box(ScrollablePane(DynamicContainer(lambda: self.current_container)), width=D(weight=1)),
                ], padding=1),
                Box(HSplit([self.footer, self.footer_hint]), height=2, padding=0, style='class:footer'),
            ]),
            floats=[
                Float(xcursor=True, ycursor=True, content=CompletionsMenu(max_height=8)),
            ],
        )
        self.application = Application(
            full_screen=True,
            layout=Layout(self.root),
            key_bindings=self.kb,
            style=PT_APP_STYLE,
            mouse_support=True,
        )

    def _initial_state(self):
        run_num = _first_available_num(self.args.recon_tag, self.args.run_tag, self.args.recon_prefix)
        return {
            'tag': self.args.recon_tag,
            'prefix': self.args.recon_prefix,
            'run_num': str(run_num),
            'config_style': 'commented',
            'workflow': 'simulation',
            'sim_model_source': 'local_pdb',
            'sim_pdb_code': SIM_DEFAULTS['make_densities']['pdb_code'],
            'sim_pdb_file': op.join(self._planned_recon_dir_for(self.args.recon_tag,
                                                                self.args.recon_prefix,
                                                                str(run_num)),
                                    'aux', '4BED.pdb'),
            'sim_detd': SIM_DEFAULTS['parameters']['detd'],
            'sim_lambda': SIM_DEFAULTS['parameters']['lambda'],
            'sim_detsize': SIM_DEFAULTS['parameters']['detsize'],
            'sim_pixsize': SIM_DEFAULTS['parameters']['pixsize'],
            'sim_stoprad': SIM_DEFAULTS['parameters']['stoprad'],
            'sim_polarization': SIM_DEFAULTS['parameters']['polarization'],
            'sim_num_data': SIM_DEFAULTS['make_data']['num_data'],
            'sim_fluence': SIM_DEFAULTS['make_data']['fluence'],
            'sim_recon_type': '3d',
            'sim_need_scaling': True,
            'sim_num_div': SIM_DEFAULTS['emc']['num_div'],
            'sim_num_rot': RECON_DEFAULTS_2D['num_rot'],
            'sim_num_modes': RECON_DEFAULTS_2D['num_modes'],
            'sim_beta_factor': SIM_DEFAULTS['emc']['beta_factor'],
            'sim_beta_jump': SIM_DEFAULTS['emc']['beta_schedule'].split()[0],
            'sim_beta_period': SIM_DEFAULTS['emc']['beta_schedule'].split()[1],
            'exp_photon_mode': 'generate_list',
            'exp_photon_file': '',
            'exp_photon_list': '',
            'exp_photon_folder': op.realpath('.'),
            'exp_photon_selected': [],
            'exp_photon_list_name': 'photons.lst',
            'exp_detector_file': '',
            'exp_output_folder': EXP_DEFAULTS['emc']['output_folder'],
            'exp_log_file': EXP_DEFAULTS['emc']['log_file'],
            'exp_recon_type': '3d',
            'exp_need_scaling': True,
            'exp_num_div': EXP_DEFAULTS['emc']['num_div'],
            'exp_num_rot': RECON_DEFAULTS_2D['num_rot'],
            'exp_num_modes': RECON_DEFAULTS_2D['num_modes'],
            'exp_beta_factor': EXP_DEFAULTS['emc']['beta_factor'],
            'exp_beta_jump': EXP_DEFAULTS['emc']['beta_schedule'].split()[0],
            'exp_beta_period': EXP_DEFAULTS['emc']['beta_schedule'].split()[1],
        }

    # Navigation, footer state, and shared widget helpers.

    def _make_key_bindings(self):
        kb = KeyBindings()

        def _focused_path_field(event):
            widget = event.app.layout.current_control
            if widget is None:
                return None
            return getattr(widget, '_path_owner', None)

        def _path_completion_active():
            return self.edit_mode and self._focused_path_widget() is not None

        @kb.add('c-c')
        @kb.add('c-q')
        def _(event):
            self.cancelled = True
            event.app.exit(result=None)

        @kb.add('tab', filter=Condition(lambda: not self.edit_mode))
        def _(event):
            event.app.layout.focus_next()

        @kb.add('tab', filter=Condition(_path_completion_active))
        def _(event):
            buffer = event.app.current_buffer
            if buffer.complete_state is None:
                buffer.start_completion(select_first=False)
            else:
                buffer.complete_next()

        @kb.add('s-tab', filter=Condition(lambda: not self.edit_mode))
        def _(event):
            event.app.layout.focus_previous()

        @kb.add('s-tab', filter=Condition(_path_completion_active))
        def _(event):
            buffer = event.app.current_buffer
            if buffer.complete_state is not None:
                buffer.complete_previous()

        @kb.add('enter')
        def _(event):
            path_field = _focused_path_field(event)
            if path_field is not None:
                self.edit_mode = not self.edit_mode
                if self.edit_mode:
                    event.app.current_buffer.start_completion(select_first=False)
                else:
                    event.app.current_buffer.cancel_completion()
                self.application.invalidate()
                return

        @kb.add('f8')
        def _(event):
            self._go_back()

        @kb.add('f9')
        def _(event):
            self._advance()

        return kb

    def _step_active(self, step_id):
        workflow = self.state['workflow']
        if step_id.startswith('sim_'):
            return workflow == 'simulation'
        if step_id.startswith('exp_'):
            return workflow == 'experimental'
        return True

    def _active_steps(self):
        return [step for step in self.STEPS if self._step_active(step)]

    def _current_index(self):
        return self._active_steps().index(self.current_step)

    def _next_step(self):
        steps = self._active_steps()
        index = self._current_index()
        if index == len(steps) - 1:
            return None
        return steps[index + 1]

    def _previous_step(self):
        steps = self._active_steps()
        index = self._current_index()
        if index == 0:
            return None
        return steps[index - 1]

    def _planned_recon_dir_for(self, tag, prefix, run_num):
        return op.join(prefix, _name_recon_dir(tag, int(run_num)))

    def _planned_recon_dir(self):
        return self._planned_recon_dir_for(self.state['tag'], self.state['prefix'], self.state['run_num'])

    def _set_status(self, message):
        self.status = message
        self.status_style = 'class:footer-info'
        self.application.invalidate()

    def _set_error_status(self, message):
        self.status = message
        self.status_style = 'class:footer-error'
        self.application.invalidate()

    def _focused_path_widget(self):
        if not hasattr(self, 'application') or self.application is None:
            return None
        control = self.application.layout.current_control
        if control is None:
            return None
        return getattr(control, '_path_owner', None)

    def _footer_mode_text(self):
        if self.edit_mode:
            return FormattedText([
                ('class:footer-mode-edit', ' EDIT MODE '),
                ('class:footer', '  '),
                ('class:footer-key', 'Tab'),
                ('class:footer', ': complete path  '),
                ('class:footer-key', 'Enter'),
                ('class:footer', ': leave edit mode'),
            ])
        if self._focused_path_widget() is not None:
            return FormattedText([
                ('class:footer-mode-path', ' PATH FIELD '),
                ('class:footer', '  '),
                ('class:footer-key', 'Enter'),
                ('class:footer', ': edit path  '),
                ('class:footer-key', 'Tab/Shift-Tab'),
                ('class:footer', ': move focus'),
            ])
        return FormattedText([
            ('class:footer-key', 'Tab/Shift-Tab'),
            ('class:footer', ': move focus  '),
            ('class:footer-key', 'F8'),
            ('class:footer', ': back  '),
            ('class:footer-key', 'F9'),
            ('class:footer', ': next'),
        ])

    def _footer_status_text(self):
        return FormattedText([
            (self.status_style, self.status),
        ])

    def _update_directory_preview(self):
        if not hasattr(self, 'widgets'):
            return
        tag = self.widgets['tag'].text.strip() or self.state['tag']
        prefix = self.widgets['prefix'].text.strip() or self.state['prefix']
        prefix = op.realpath(prefix)
        run_value = self.widgets['run_num'].text.strip() or self.state['run_num']
        try:
            current_run = int(run_value)
        except ValueError:
            current_run = int(self.state['run_num'])
        try:
            old_auto = _first_available_num(self.state['tag'], int(self.state['run_num']), self.state['prefix'])
        except ValueError:
            old_auto = int(self.state['run_num'])
        if run_value == str(self.state['run_num']) or current_run == old_auto:
            new_auto = _first_available_num(tag, max(1, current_run), prefix)
            self.widgets['run_num'].text = str(new_auto)
            current_run = new_auto
        self.directory_preview.text = 'Current proposed directory: %s' % self._planned_recon_dir_for(tag, prefix, current_run)

    def _summary_entries(self):
        entries = []

        def add(key, label, value, step_id):
            if value in (None, ''):
                return
            entries.append((key, '%s: %s' % (label, value), step_id))

        add('directory', 'Reconstruction directory', self._planned_recon_dir(), 'directory')
        add('config_style', 'Config style', self.state['config_style'], 'config_style')
        add('workflow', 'Workflow', self.state['workflow'], 'workflow')
        if self.state['workflow'] == 'simulation':
            add('sim_model_source', 'Simulation model source', self.state['sim_model_source'], 'sim_model_source')
            if self.state['sim_model_source'] == 'pdb_code':
                add('sim_pdb_code', 'PDB code', self.state['sim_pdb_code'], 'sim_model_value')
            else:
                add('sim_pdb_file', 'Local PDB file', self.state['sim_pdb_file'], 'sim_model_value')
            add('sim_detd', 'Detector distance', self.state['sim_detd'], 'sim_geometry')
            add('sim_lambda', 'Photon wavelength', self.state['sim_lambda'], 'sim_geometry')
            add('sim_detsize', 'Detector size', self.state['sim_detsize'], 'sim_geometry')
            add('sim_pixsize', 'Pixel size', self.state['sim_pixsize'], 'sim_geometry')
            add('sim_stoprad', 'Beamstop radius', self.state['sim_stoprad'], 'sim_geometry')
            add('sim_polarization', 'Polarization', self.state['sim_polarization'], 'sim_geometry')
            add('sim_num_data', 'Number of patterns', self.state['sim_num_data'], 'sim_data')
            add('sim_fluence', 'Fluence', self.state['sim_fluence'], 'sim_data')
            add('sim_recon_type', 'Recon type', self.state['sim_recon_type'], 'sim_recon')
            add('sim_need_scaling', 'Need scaling', '1' if self.state['sim_need_scaling'] else '0', 'sim_recon')
            if self.state['sim_recon_type'] == '3d':
                add('sim_num_div', 'num_div', self.state['sim_num_div'], 'sim_recon')
            else:
                add('sim_num_rot', 'num_rot', self.state['sim_num_rot'], 'sim_recon')
                add('sim_num_modes', 'num_modes', self.state['sim_num_modes'], 'sim_recon')
            add('sim_beta_factor', 'beta_factor', self.state['sim_beta_factor'], 'sim_beta')
            add('sim_beta_schedule', 'beta_schedule',
                '%s %s' % (self.state['sim_beta_jump'], self.state['sim_beta_period']), 'sim_beta')
        else:
            add('exp_photon_mode', 'Photon input mode', self.state['exp_photon_mode'], 'exp_photon_mode')
            if self.state['exp_photon_mode'] == 'single_file':
                add('exp_photon_file', 'Photon file', self.state['exp_photon_file'], 'exp_photon_input')
            elif self.state['exp_photon_mode'] == 'list_file':
                add('exp_photon_list', 'Photon list', self.state['exp_photon_list'], 'exp_photon_input')
            else:
                add('exp_photon_folder', 'Photon folder', self.state['exp_photon_folder'], 'exp_photon_input')
                if self.state['exp_photon_selected']:
                    add('exp_photon_selected', 'Selected photon files',
                        ', '.join(self.state['exp_photon_selected']), 'exp_photon_input')
                add('exp_photon_list_name', 'Photon list filename',
                    self.state['exp_photon_list_name'], 'exp_photon_input')
            add('exp_detector_file', 'Detector file', self.state['exp_detector_file'], 'exp_detector')
            add('exp_output_folder', 'Output folder', self.state['exp_output_folder'], 'exp_detector')
            add('exp_log_file', 'Log file', self.state['exp_log_file'], 'exp_detector')
            add('exp_recon_type', 'Recon type', self.state['exp_recon_type'], 'exp_recon')
            add('exp_need_scaling', 'Need scaling', '1' if self.state['exp_need_scaling'] else '0', 'exp_recon')
            if self.state['exp_recon_type'] == '3d':
                add('exp_num_div', 'num_div', self.state['exp_num_div'], 'exp_recon')
            else:
                add('exp_num_rot', 'num_rot', self.state['exp_num_rot'], 'exp_recon')
                add('exp_num_modes', 'num_modes', self.state['exp_num_modes'], 'exp_recon')
            add('exp_beta_factor', 'beta_factor', self.state['exp_beta_factor'], 'exp_beta')
            add('exp_beta_schedule', 'beta_schedule',
                '%s %s' % (self.state['exp_beta_jump'], self.state['exp_beta_period']), 'exp_beta')
        return entries

    def _build_summary_container(self):
        entries = self._summary_entries()
        if not entries:
            body = self._make_label('No selections yet')
            return Frame(Box(body, padding=1), title='Selections', width=D(preferred=44))
        labels = []
        for _, text, step_id in entries:
            style = 'class:focused' if step_id == self.current_step else 'class:label'
            labels.append(Label(text=text, style=style))
        return Frame(Box(HSplit(labels, padding=0), padding=1), title='Selections', width=D(preferred=44))

    # Step builders and validators, in wizard order.

    def _make_text(self, value='', path=False, directories=False):
        kwargs = {'text': str(value), 'multiline': False, 'style': 'class:input-field'}
        if path:
            kwargs['completer'] = PathCompleter(only_directories=directories, expanduser=True)
            kwargs['complete_while_typing'] = True
        field = TextArea(**kwargs)
        if path:
            field.buffer.read_only = Condition(
                lambda field=field: not (self.edit_mode and self.application.layout.has_focus(field))
            )
        field.control._path_owner = field if path else None
        field._is_path_field = path
        return field

    def _make_label(self, text):
        return Label(text=text, style='class:label')

    def _make_radio(self, values, current_value):
        radio = RadioList(values)
        radio.current_value = current_value
        return radio

    def _buttons(self, extra=None):
        buttons = []
        if extra is not None:
            buttons.extend(extra)
        next_text = 'Finish (F9)' if self._next_step() is None else 'Next (F9)'
        buttons.append(Button(next_text, handler=self._advance, left_symbol='', right_symbol=''))
        if self._previous_step() is not None:
            buttons.append(Button('Back (F8)', handler=self._go_back, left_symbol='', right_symbol=''))
        return HSplit(buttons, padding=1)

    def _go_back(self):
        prev_step = self._previous_step()
        if prev_step is not None:
            self.edit_mode = False
            self.current_step = prev_step
            self._set_status('Moved back to %s.' % self.STEPS[self.current_step])
            self._rebuild_ui()

    def _clear_workflow_state(self, workflow):
        defaults = self._initial_state()
        if workflow == 'simulation':
            for key in list(self.state):
                if key.startswith('exp_'):
                    self.state[key] = defaults[key]
        else:
            for key in list(self.state):
                if key.startswith('sim_'):
                    self.state[key] = defaults[key]
            self.state['sim_pdb_file'] = op.join(self._planned_recon_dir(), 'aux', '4BED.pdb')

    # Summary pane.

    def _build_directory_step(self):
        tag = self._make_text(self.state['tag'])
        prefix = self._make_text(self.state['prefix'], path=True, directories=True)
        run_num = self._make_text(self.state['run_num'])
        preview = self._make_label('Current proposed directory: %s' % self._planned_recon_dir())
        tag.buffer.on_text_changed += lambda _: self._update_directory_preview()
        prefix.buffer.on_text_changed += lambda _: self._update_directory_preview()
        run_num.buffer.on_text_changed += lambda _: self._update_directory_preview()
        self.directory_preview = preview
        self.current_focus = tag
        self.widgets = {'tag': tag, 'prefix': prefix, 'run_num': run_num}
        return Frame(Box(HSplit([
            self._make_label('Set the reconstruction tag, parent directory, and run number.'),
            self._make_label('Reconstruction tag'), tag,
            self._make_label('Parent directory'), prefix,
            self._make_label('Run number'), run_num,
            preview,
            self._buttons(),
        ]), padding=1), title=self.STEPS[self.current_step])

    def _validate_directory_step(self):
        tag = self.widgets['tag'].text.strip()
        prefix = op.realpath(self.widgets['prefix'].text.strip())
        run_num = self.widgets['run_num'].text.strip()
        if not tag:
            self._set_error_status('Reconstruction tag is required.')
            return False
        if not op.isdir(prefix):
            self._set_error_status('Parent directory does not exist: %s' % prefix)
            return False
        try:
            run_val = int(run_num)
        except ValueError:
            self._set_error_status('Run number must be an integer.')
            return False
        if run_val < 1:
            self._set_error_status('Run number must be at least 1.')
            return False
        recon_dir = self._planned_recon_dir_for(tag, prefix, run_num)
        if op.exists(recon_dir):
            self._set_error_status('That reconstruction directory already exists: %s' % recon_dir)
            return False
        self.state['tag'] = tag
        self.state['prefix'] = prefix
        self.state['run_num'] = str(run_val)
        self.state['sim_pdb_file'] = op.join(recon_dir, 'aux', '4BED.pdb')
        return True

    def _build_config_style_step(self):
        radio = self._make_radio([('commented', 'Keep helpful comments'),
                                  ('clean', 'Write a clean config')],
                                 self.state['config_style'])
        self.current_focus = radio
        self.widgets = {'radio': radio}
        return Frame(Box(HSplit([
            Label(text='Choose how much commentary to include in the generated config file.'),
            radio,
            self._buttons(),
        ]), padding=1), title=self.STEPS[self.current_step])

    def _validate_config_style_step(self):
        self.state['config_style'] = self.widgets['radio'].current_value
        return True

    def _build_workflow_step(self):
        radio = self._make_radio([('simulation', 'Simulation'), ('experimental', 'Experimental')],
                                 self.state['workflow'])
        self.current_focus = radio
        self.widgets = {'radio': radio}
        return Frame(Box(HSplit([
            Label(text='Choose whether to generate a simulation config or an experimental EMC config.'),
            radio,
            self._buttons(),
        ]), padding=1), title=self.STEPS[self.current_step])

    def _validate_workflow_step(self):
        workflow = self.widgets['radio'].current_value
        changed = workflow != self.state['workflow']
        self.state['workflow'] = workflow
        if changed:
            self._clear_workflow_state(workflow)
        return True

    def _build_sim_model_source_step(self):
        radio = self._make_radio([('pdb_code', 'Fetch by PDB code'), ('local_pdb', 'Use a local PDB file')],
                                 self.state['sim_model_source'])
        self.current_focus = radio
        self.widgets = {'radio': radio}
        return Frame(Box(HSplit([
            Label(text='Select whether the simulation should start from a PDB code or a local PDB file.'),
            radio,
            self._buttons(),
        ]), padding=1), title=self.STEPS[self.current_step])

    def _validate_sim_model_source_step(self):
        self.state['sim_model_source'] = self.widgets['radio'].current_value
        return True

    def _build_sim_model_value_step(self):
        if self.state['sim_model_source'] == 'pdb_code':
            field = self._make_text(self.state['sim_pdb_code'])
            label = 'PDB code'
            self.widgets = {'field': field}
            self.current_focus = field
        else:
            field = self._make_text(self.state['sim_pdb_file'], path=True)
            label = 'Local PDB file path'
            self.widgets = {'field': field}
            self.current_focus = field
        return Frame(Box(HSplit([
            Label(text='Provide the structure input for the simulation.'),
            Label(text=label),
            field,
            self._buttons(),
        ]), padding=1), title=self.STEPS[self.current_step])

    def _validate_sim_model_value_step(self):
        value = self.widgets['field'].text.strip()
        if self.state['sim_model_source'] == 'pdb_code':
            if not value:
                self._set_error_status('PDB code is required.')
                return False
            self.state['sim_pdb_code'] = value
            return True
        recon_local = op.join(self._planned_recon_dir(), 'aux', '4BED.pdb')
        install_local = op.join(self.parent_dir, 'aux', '4BED.pdb')
        if not value:
            self._set_error_status('Local PDB file path is required.')
            return False
        if op.realpath(value) == op.realpath(recon_local) and op.exists(install_local):
            self.state['sim_pdb_file'] = recon_local
            return True
        if not op.exists(value):
            self._set_error_status('PDB file does not exist: %s' % value)
            return False
        self.state['sim_pdb_file'] = op.realpath(value)
        return True

    def _build_sim_geometry_step(self):
        widgets = {
            'detd': self._make_text(self.state['sim_detd']),
            'lambda': self._make_text(self.state['sim_lambda']),
            'detsize': self._make_text(self.state['sim_detsize']),
            'pixsize': self._make_text(self.state['sim_pixsize']),
            'stoprad': self._make_text(self.state['sim_stoprad']),
            'polarization': self._make_radio([('x', 'x'), ('y', 'y'), ('none', 'none')],
                                             self.state['sim_polarization']),
        }
        self.widgets = widgets
        self.current_focus = widgets['detd']
        return Frame(Box(HSplit([
            Label(text='Enter detector geometry values for the simulation.'),
            Label(text='Detector distance (mm)'), widgets['detd'],
            Label(text='Photon wavelength (Angstrom)'), widgets['lambda'],
            Label(text='Detector size (one value or X Y)'), widgets['detsize'],
            Label(text='Pixel size (mm)'), widgets['pixsize'],
            Label(text='Beamstop radius (pixels)'), widgets['stoprad'],
            Label(text='Polarization correction'), widgets['polarization'],
            self._buttons(),
        ]), padding=1), title=self.STEPS[self.current_step])

    def _validate_float_text(self, value, label, minimum=0):
        try:
            parsed = float(value)
        except ValueError:
            self._set_error_status('%s must be numeric.' % label)
            return None
        if parsed < minimum:
            self._set_error_status('%s must be at least %s.' % (label, minimum))
            return None
        return value

    def _validate_int_text(self, value, label, minimum=1):
        try:
            parsed = int(value)
        except ValueError:
            self._set_error_status('%s must be an integer.' % label)
            return None
        if parsed < minimum:
            self._set_error_status('%s must be at least %d.' % (label, minimum))
            return None
        return str(parsed)

    def _validate_detsize_text(self, value):
        parts = value.split()
        if len(parts) not in (1, 2):
            self._set_error_status('Detector size must be one integer or two integers.')
            return None
        try:
            parsed = [int(part) for part in parts]
        except ValueError:
            self._set_error_status('Detector size must contain integers only.')
            return None
        if min(parsed) <= 0:
            self._set_error_status('Detector size values must be positive.')
            return None
        return ' '.join(str(item) for item in parsed)

    def _validate_sim_geometry_step(self):
        detd = self._validate_float_text(self.widgets['detd'].text.strip(), 'Detector distance')
        wavelength = self._validate_float_text(self.widgets['lambda'].text.strip(), 'Photon wavelength')
        detsize = self._validate_detsize_text(self.widgets['detsize'].text.strip())
        pixsize = self._validate_float_text(self.widgets['pixsize'].text.strip(), 'Pixel size')
        stoprad = self._validate_float_text(self.widgets['stoprad'].text.strip(), 'Beamstop radius')
        if None in (detd, wavelength, detsize, pixsize, stoprad):
            return False
        self.state['sim_detd'] = detd
        self.state['sim_lambda'] = wavelength
        self.state['sim_detsize'] = detsize
        self.state['sim_pixsize'] = pixsize
        self.state['sim_stoprad'] = stoprad
        self.state['sim_polarization'] = self.widgets['polarization'].current_value
        return True

    def _build_sim_data_step(self):
        num_data = self._make_text(self.state['sim_num_data'])
        fluence = self._make_text(self.state['sim_fluence'])
        self.widgets = {'num_data': num_data, 'fluence': fluence}
        self.current_focus = num_data
        return Frame(Box(HSplit([
            Label(text='Set the number of simulated diffraction patterns and fluence.'),
            Label(text='Number of diffraction patterns'), num_data,
            Label(text='Incident fluence (photons/um^2)'), fluence,
            self._buttons(),
        ]), padding=1), title=self.STEPS[self.current_step])

    def _validate_sim_data_step(self):
        num_data = self._validate_int_text(self.widgets['num_data'].text.strip(), 'Number of patterns')
        fluence = self._validate_float_text(self.widgets['fluence'].text.strip(), 'Fluence')
        if None in (num_data, fluence):
            return False
        self.state['sim_num_data'] = num_data
        self.state['sim_fluence'] = fluence
        return True

    def _switch_recon_type(self, prefix, recon_type):
        self.state['%s_recon_type' % prefix] = recon_type
        if 'scaling' in self.widgets:
            self.state['%s_need_scaling' % prefix] = self.widgets['scaling'].current_value
        if 'num_div' in self.widgets:
            self.state['%s_num_div' % prefix] = self.widgets['num_div'].text.strip() or self.state['%s_num_div' % prefix]
        if 'num_rot' in self.widgets:
            self.state['%s_num_rot' % prefix] = self.widgets['num_rot'].text.strip() or self.state['%s_num_rot' % prefix]
        if 'num_modes' in self.widgets:
            self.state['%s_num_modes' % prefix] = self.widgets['num_modes'].text.strip() or self.state['%s_num_modes' % prefix]
        self._set_status('Reconstruction type set to %s.' % recon_type)
        self._rebuild_ui()

    def _recon_button_label(self, current_recon, candidate):
        if current_recon == candidate:
            return '%s selected' % candidate.upper()
        return candidate.upper()

    def _build_recon_settings_step(self, prefix):
        scaling = self._make_radio([(True, 'Enable fluence scaling'), (False, 'Disable fluence scaling')],
                                   self.state['%s_need_scaling' % prefix])
        current_recon = self.state['%s_recon_type' % prefix]
        button_3d = Button(self._recon_button_label(current_recon, '3d'),
                           handler=lambda: self._switch_recon_type(prefix, '3d'),
                           left_symbol='', right_symbol='')
        button_2d = Button(self._recon_button_label(current_recon, '2d'),
                           handler=lambda: self._switch_recon_type(prefix, '2d'),
                           left_symbol='', right_symbol='')
        widgets = {'scaling': scaling, 'button_3d': button_3d, 'button_2d': button_2d}
        parts = [
            self._make_label('Set reconstruction type and scaling options.'),
            self._make_label('Reconstruction type'),
            VSplit([button_3d, button_2d], padding=1),
            self._make_label('Current reconstruction type: %s' % current_recon.upper()),
            self._make_label('Fluence scaling'), scaling,
        ]
        if current_recon == '3d':
            num_div = self._make_text(self.state['%s_num_div' % prefix])
            widgets['num_div'] = num_div
            parts.extend([self._make_label('Quaternion sampling num_div'), num_div])
        else:
            num_rot = self._make_text(self.state['%s_num_rot' % prefix])
            num_modes = self._make_text(self.state['%s_num_modes' % prefix])
            widgets['num_rot'] = num_rot
            widgets['num_modes'] = num_modes
            parts.extend([
                self._make_label('Number of in-plane rotations'), num_rot,
                self._make_label('Number of modes'), num_modes,
            ])
        parts.append(self._buttons())
        self.widgets = widgets
        self.current_focus = button_3d
        return Frame(Box(HSplit(parts), padding=1), title=self.STEPS[self.current_step])

    def _validate_recon_settings_step(self, prefix):
        new_type = self.state['%s_recon_type' % prefix]
        self.state['%s_need_scaling' % prefix] = self.widgets['scaling'].current_value
        if new_type == '3d':
            num_div = self._validate_int_text(self.widgets['num_div'].text.strip(), 'num_div')
            if num_div is None:
                return False
            self.state['%s_num_div' % prefix] = num_div
        else:
            num_rot = self._validate_int_text(self.widgets['num_rot'].text.strip(), 'num_rot')
            num_modes = self._validate_int_text(self.widgets['num_modes'].text.strip(), 'num_modes')
            if None in (num_rot, num_modes):
                return False
            self.state['%s_num_rot' % prefix] = num_rot
            self.state['%s_num_modes' % prefix] = num_modes
        return True

    def _build_sim_recon_step(self):
        return self._build_recon_settings_step('sim')

    def _validate_sim_recon_step(self):
        return self._validate_recon_settings_step('sim')

    def _build_beta_step(self, prefix):
        factor = self._make_text(self.state['%s_beta_factor' % prefix])
        jump = self._make_text(self.state['%s_beta_jump' % prefix])
        period = self._make_text(self.state['%s_beta_period' % prefix])
        self.widgets = {'factor': factor, 'jump': jump, 'period': period}
        self.current_focus = factor
        return Frame(Box(HSplit([
            Label(text='beta_start[d] is computed per frame. The iteration factor is'),
            Label(text='beta_factor * beta_schedule[0]**((i-1)//beta_schedule[1]).'),
            Label(text='Initial beta_factor'), factor,
            Label(text='Multiply beta_factor by how much at each change?'), jump,
            Label(text='Change beta_factor every how many iterations?'), period,
            self._buttons(),
        ]), padding=1), title=self.STEPS[self.current_step])

    def _validate_beta_step(self, prefix):
        factor = self._validate_float_text(self.widgets['factor'].text.strip(), 'beta_factor')
        jump = self._validate_float_text(self.widgets['jump'].text.strip(), 'beta_jump')
        period = self._validate_int_text(self.widgets['period'].text.strip(), 'beta_period')
        if None in (factor, jump, period):
            return False
        self.state['%s_beta_factor' % prefix] = factor
        self.state['%s_beta_jump' % prefix] = jump
        self.state['%s_beta_period' % prefix] = period
        return True

    def _build_sim_beta_step(self):
        return self._build_beta_step('sim')

    def _validate_sim_beta_step(self):
        return self._validate_beta_step('sim')

    def _build_exp_photon_mode_step(self):
        radio = self._make_radio([
            ('single_file', 'Use a single existing photon file'),
            ('list_file', 'Use an existing photon list file'),
            ('generate_list', 'Create a new photon list from a folder'),
        ], self.state['exp_photon_mode'])
        self.widgets = {'radio': radio}
        self.current_focus = radio
        return Frame(Box(HSplit([
            Label(text='Choose how photon inputs should be configured.'),
            radio,
            self._buttons(),
        ]), padding=1), title=self.STEPS[self.current_step])

    def _validate_exp_photon_mode_step(self):
        self.state['exp_photon_mode'] = self.widgets['radio'].current_value
        return True

    def _generated_file_values(self, folder):
        if not op.isdir(folder):
            return [(EMPTY_CHECKBOX_VALUE, 'No files found')]
        entries = [name for name in sorted(os.listdir(folder)) if op.isfile(op.join(folder, name))]
        if not entries:
            return [(EMPTY_CHECKBOX_VALUE, 'No files found')]
        return [(name, name) for name in entries]

    def _sync_generated_files(self, announce=False):
        if not hasattr(self, 'widgets') or 'folder' not in self.widgets or 'checklist' not in self.widgets:
            return
        folder = self.widgets['folder'].text.strip()
        checklist = self.widgets['checklist']
        values = self._generated_file_values(folder)
        valid_names = set(name for name, _ in values if name != EMPTY_CHECKBOX_VALUE)
        selected = [name for name in checklist.current_values if name in valid_names]
        checklist.values = values
        checklist.current_values = selected
        self.state['exp_photon_folder'] = folder
        self.state['exp_photon_selected'] = selected
        if announce:
            self._set_status('Photon file list refreshed for %s.' % folder)
            self.application.invalidate()

    def _refresh_generated_files(self):
        self._sync_generated_files(announce=True)

    def _build_exp_photon_input_step(self):
        mode = self.state['exp_photon_mode']
        parts = [Label(text='Provide the photon input details for the chosen mode.')]
        widgets = {}
        extra = []
        if mode == 'single_file':
            field = self._make_text(self.state['exp_photon_file'], path=True)
            widgets['field'] = field
            parts.extend([Label(text='Photon file path'), field])
            self.current_focus = field
        elif mode == 'list_file':
            field = self._make_text(self.state['exp_photon_list'], path=True)
            widgets['field'] = field
            parts.extend([Label(text='Photon list file path'), field])
            self.current_focus = field
        else:
            folder = self._make_text(self.state['exp_photon_folder'], path=True, directories=True)
            values = self._generated_file_values(folder.text.strip())
            checklist = CheckboxList(values)
            checklist.current_values = list(self.state['exp_photon_selected'])
            list_name = self._make_text(self.state['exp_photon_list_name'])
            folder.buffer.on_text_changed += lambda _: self._sync_generated_files()
            widgets['folder'] = folder
            widgets['checklist'] = checklist
            widgets['list_name'] = list_name
            self.current_focus = folder
            extra.append(Button('Refresh Files', handler=self._refresh_generated_files,
                                left_symbol='', right_symbol=''))
            parts.extend([
                Label(text='Folder containing photon files'), folder,
                Label(text='Select files to include'), checklist,
                Label(text='Photon list filename to write later'), list_name,
            ])
        parts.append(self._buttons(extra=extra))
        self.widgets = widgets
        return Frame(Box(HSplit(parts), padding=1), title=self.STEPS[self.current_step])

    def _validate_exp_photon_input_step(self):
        mode = self.state['exp_photon_mode']
        if mode == 'single_file':
            value = self.widgets['field'].text.strip()
            if not op.exists(value):
                self._set_error_status('Photon file does not exist: %s' % value)
                return False
            self.state['exp_photon_file'] = op.realpath(value)
            self.state['exp_photon_list'] = ''
            return True
        if mode == 'list_file':
            value = self.widgets['field'].text.strip()
            if not op.exists(value):
                self._set_error_status('Photon list file does not exist: %s' % value)
                return False
            self.state['exp_photon_list'] = op.realpath(value)
            self.state['exp_photon_file'] = ''
            return True
        folder = op.realpath(self.widgets['folder'].text.strip())
        if not op.isdir(folder):
            self._set_error_status('Photon folder does not exist: %s' % folder)
            return False
        selected = list(self.widgets['checklist'].current_values)
        if not selected:
            self._set_error_status('Select at least one photon file.')
            return False
        list_name = self.widgets['list_name'].text.strip()
        if not list_name:
            self._set_error_status('Photon list filename is required.')
            return False
        self.state['exp_photon_folder'] = folder
        self.state['exp_photon_selected'] = selected
        self.state['exp_photon_list_name'] = list_name
        self.state['exp_photon_file'] = ''
        self.state['exp_photon_list'] = ''
        return True

    def _build_exp_detector_step(self):
        detector = self._make_text(self.state['exp_detector_file'], path=True)
        output_folder = self._make_text(self.state['exp_output_folder'])
        log_file = self._make_text(self.state['exp_log_file'])
        self.widgets = {'detector': detector, 'output_folder': output_folder, 'log_file': log_file}
        self.current_focus = detector
        return Frame(Box(HSplit([
            Label(text='Provide detector and output paths for the experimental setup.'),
            Label(text='Detector file path'), detector,
            Label(text='Output folder'), output_folder,
            Label(text='Log file'), log_file,
            self._buttons(),
        ]), padding=1), title=self.STEPS[self.current_step])

    def _validate_exp_detector_step(self):
        detector = self.widgets['detector'].text.strip()
        if not op.exists(detector):
            self._set_error_status('Detector file does not exist: %s' % detector)
            return False
        output_folder = self.widgets['output_folder'].text.strip()
        log_file = self.widgets['log_file'].text.strip()
        if not output_folder or not log_file:
            self._set_error_status('Output folder and log file are required.')
            return False
        self.state['exp_detector_file'] = op.realpath(detector)
        self.state['exp_output_folder'] = output_folder
        self.state['exp_log_file'] = log_file
        return True

    def _build_exp_recon_step(self):
        return self._build_recon_settings_step('exp')

    def _validate_exp_recon_step(self):
        return self._validate_recon_settings_step('exp')

    def _build_exp_beta_step(self):
        return self._build_beta_step('exp')

    def _validate_exp_beta_step(self):
        return self._validate_beta_step('exp')

    def _build_current_container(self):
        return getattr(self, '_build_%s_step' % self.current_step)()

    def _validate_current_step(self):
        return getattr(self, '_validate_%s_step' % self.current_step)()

    def _rebuild_ui(self):
        if self.edit_mode and (self.current_focus is None or not getattr(self.current_focus, '_is_path_field', False)):
            self.edit_mode = False
        self.summary_container = self._build_summary_container()
        self.current_container = self._build_current_container()
        self.application.invalidate()
        if self.current_focus is not None:
            self.application.layout.focus(self.current_focus)

    def _advance(self):
        result = self._validate_current_step()
        if result is False:
            self._rebuild_ui()
            return
        next_step = self._next_step()
        if next_step is None:
            self.application.exit(result=self.state)
            return
        self.edit_mode = False
        self.current_step = next_step
        self._rebuild_ui()

    # Running the wizard and translating state into config output.

    def run(self):
        self._rebuild_ui()
        return self.application.run()

    def _simulation_config(self):
        config = {section: values.copy() for section, values in SIM_DEFAULTS.items()}
        if self.state['sim_model_source'] == 'pdb_code':
            config['make_densities']['pdb_code'] = self.state['sim_pdb_code']
        else:
            config['make_densities'].pop('pdb_code', None)
            if op.realpath(self.state['sim_pdb_file']) == op.realpath(op.join(self.parent_dir, 'aux', '4BED.pdb')):
                config['make_densities']['in_pdb_file'] = 'aux/4BED.pdb'
            else:
                config['make_densities']['in_pdb_file'] = self.state['sim_pdb_file']
        config['parameters']['detd'] = self.state['sim_detd']
        config['parameters']['lambda'] = self.state['sim_lambda']
        config['parameters']['detsize'] = self.state['sim_detsize']
        config['parameters']['pixsize'] = self.state['sim_pixsize']
        config['parameters']['stoprad'] = self.state['sim_stoprad']
        config['parameters']['polarization'] = self.state['sim_polarization']
        config['make_data']['num_data'] = self.state['sim_num_data']
        config['make_data']['fluence'] = self.state['sim_fluence']
        emc = config['emc']
        emc['recon_type'] = self.state['sim_recon_type']
        emc['need_scaling'] = '1' if self.state['sim_need_scaling'] else '0'
        emc.pop('num_div', None)
        emc.pop('num_rot', None)
        emc.pop('num_modes', None)
        if self.state['sim_recon_type'] == '3d':
            emc['num_div'] = self.state['sim_num_div']
        else:
            emc['num_rot'] = self.state['sim_num_rot']
            emc['num_modes'] = self.state['sim_num_modes']
        emc['beta_factor'] = self.state['sim_beta_factor']
        emc['beta_schedule'] = '%s %s' % (self.state['sim_beta_jump'], self.state['sim_beta_period'])
        return config

    def _experimental_config(self):
        config = {section: values.copy() for section, values in EXP_DEFAULTS.items()}
        emc = config['emc']
        mode = self.state['exp_photon_mode']
        if mode == 'single_file':
            emc['in_photons_file'] = self.state['exp_photon_file']
        else:
            emc['in_photons_list'] = self.state['exp_photon_list'] if mode == 'list_file' else self.state['exp_photon_list_name']
        emc['in_detector_file'] = self.state['exp_detector_file']
        emc['recon_type'] = self.state['exp_recon_type']
        emc['output_folder'] = self.state['exp_output_folder']
        emc['log_file'] = self.state['exp_log_file']
        emc['need_scaling'] = '1' if self.state['exp_need_scaling'] else '0'
        emc.pop('num_div', None)
        emc.pop('num_rot', None)
        emc.pop('num_modes', None)
        if self.state['exp_recon_type'] == '3d':
            emc['num_div'] = self.state['exp_num_div']
        else:
            emc['num_rot'] = self.state['exp_num_rot']
            emc['num_modes'] = self.state['exp_num_modes']
        emc['beta_factor'] = self.state['exp_beta_factor']
        emc['beta_schedule'] = '%s %s' % (self.state['exp_beta_jump'], self.state['exp_beta_period'])
        config['emc'] = _ordered_experimental_emc(emc)
        return config

    def write_generated_files(self):
        recon_dir = _create_recon_dir(tag=self.state['tag'],
                                      num=int(self.state['run_num']),
                                      prefix=self.state['prefix'])
        _setup_aux_dir(recon_dir, self.parent_dir, copy_aux=self.args.copy_aux)
        if self.state['workflow'] == 'experimental' and self.state['exp_photon_mode'] == 'generate_list':
            list_path = op.join(recon_dir, self.state['exp_photon_list_name'])
            with open(list_path, 'w') as fptr:
                for name in self.state['exp_photon_selected']:
                    fptr.write(op.join(self.state['exp_photon_folder'], name) + '\n')
        config = self._simulation_config() if self.state['workflow'] == 'simulation' else self._experimental_config()
        keep_comments = self.state['config_style'] == 'commented'
        workflow_name = 'simulation' if self.state['workflow'] == 'simulation' else 'experimental'
        _write_generated_config(recon_dir, _render_config(config, keep_comments, workflow_name))
        return recon_dir

# Config rendering and command-line entry points.

def _render_config(config, keep_comments, workflow):
    lines = []
    if keep_comments:
        lines.extend([
            '# Generated by dragonfly.init',
            '# Edit this file later if you want to refine the setup.',
            '# For details about options, see:',
            '#     https://github.com/duaneloh/Dragonfly/wiki/Configuring-your-experiment',
            '',
        ])
        if workflow == 'experimental':
            lines.extend([
                '# Experimental workflows often only need the [emc] section here.',
                '# Add [parameters] and [make_detector] later if you need detector generation.',
                '',
            ])

    for section, values in config.items():
        if keep_comments:
            if section == 'parameters':
                lines.append('# Detector geometry in mm / Angstrom / pixels')
            elif section == 'make_densities':
                lines.append('# Density map setup for simulated data')
            elif section == 'make_data':
                lines.append('# Simulated photon data generation')
            elif section == 'emc':
                lines.extend([
                    '# beta_start[d] is computed per frame. For iteration i,',
                    '# factor = beta_factor * beta_schedule[0]**((i-1)//beta_schedule[1])',
                    '# and beta[d] = beta_start[d] * factor',
                ])
        lines.append('[%s]' % section)
        for key, value in values.items():
            lines.append('%s = %s' % (key, value))
        lines.append('')
    return '\n'.join(lines).rstrip() + '\n'

def _write_generated_config(recon_dir, config_text):
    with open(op.join(recon_dir, 'config.ini'), 'w') as fptr:
        fptr.write(config_text)

def _run_legacy_setup(args, parent_dir):
    new_recon_dir = _create_recon_dir(tag=args.recon_tag, num=args.run_tag,
                                      prefix=args.recon_prefix)
    print(80 * '=')
    print('Initializing new directory and creating soft links to useful utilities.')
    print("Type 'dragonfly.init -h' for options")
    print('See https://dragonfly-spi.readthedocs.io/en/latest/user-guides/faq.html for troubleshooting tips.')
    print(80 * '=')
    if args.recon_prefix != './':
        print('Created new directory:', new_recon_dir)
    _setup_aux_dir(new_recon_dir, parent_dir, copy_aux=args.copy_aux)
    _copy_default_config(new_recon_dir, parent_dir)

def _run_interactive_setup(args, parent_dir):
    if not HAS_INTERACTIVE_TUI:
        print('prompt_toolkit/rich unavailable; falling back to --defaults behavior.')
        _run_legacy_setup(args, parent_dir)
        return

    wizard = FullScreenWizard(args, parent_dir)
    result = wizard.run()
    if wizard.cancelled or result is None:
        raise KeyboardInterrupt
    recon_dir = wizard.write_generated_files()
    console = Console()
    console.rule('[bold cyan]Setup Complete[/bold cyan]')
    console.print('Created new directory: ', end='')
    console.print(recon_dir, style='green', markup=False)
    console.print('Config file: ', end='')
    console.print(op.join(recon_dir, 'config.ini'), style='green', markup=False)
    console.print('Next step: ', end='')
    console.print('cd %s' % recon_dir, style='green', markup=False)

def main():
    '''Parse command line arguments and create a new reconstruction directory.'''
    parser = argparse.ArgumentParser(
        'Creates a new reconstruction instance based on the packaged template',
    )
    parser.add_argument('-t', '--recon_file_tag', dest='recon_tag', default='recon',
                        help='Prefix the reconstruction folders with your specified tag.')
    parser.add_argument('-r', '--recon_run_num', dest='run_tag', type=int, default=1,
                        help='Give your reconstruction a specific number if it does not already exist.')
    parser.add_argument('-p', '--recon_prefix', dest='recon_prefix', default='./',
                        help='Path to the folder containing the reconstruction folder.')
    parser.add_argument('--defaults', '--legacy', dest='legacy_mode', action='store_true',
                        help='Keep the old non-interactive behavior and copy the packaged config.ini.')
    parser.add_argument('--copy-aux', dest='copy_aux', action='store_true',
                        help='Copy aux/ into the reconstruction directory instead of symlinking it.')
    args = parser.parse_args()
    args.recon_prefix = op.realpath(args.recon_prefix)

    parent_dir = op.realpath(op.dirname(op.dirname(op.realpath(__file__))))
    interactive = sys.stdin.isatty() and sys.stdout.isatty() and not args.legacy_mode
    if interactive:
        _run_interactive_setup(args, parent_dir)
    else:
        _run_legacy_setup(args, parent_dir)

if __name__ == '__main__':
    main()
