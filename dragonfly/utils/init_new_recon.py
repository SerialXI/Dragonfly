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
    from prompt_toolkit import prompt as pt_prompt
    from prompt_toolkit.styles import Style as PTStyle
    from prompt_toolkit.completion import PathCompleter
    from prompt_toolkit.shortcuts import checkboxlist_dialog, radiolist_dialog, yes_no_dialog
    HAS_PROMPT_TOOLKIT = True
except ImportError:
    PTStyle = None
    PathCompleter = None
    checkboxlist_dialog = None
    radiolist_dialog = None
    yes_no_dialog = None
    pt_prompt = None
    HAS_PROMPT_TOOLKIT = False

try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    HAS_RICH = True
except ImportError:
    Console = None
    Panel = None
    Table = None
    HAS_RICH = False

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

COLOR = {
    'title': '\033[95m',
    'section': '\033[94m',
    'prompt': '\033[96m',
    'value': '\033[92m',
    'warn': '\033[93m',
    'error': '\033[91m',
    'reset': '\033[0m',
    'bold': '\033[1m',
}

CONSOLE = Console() if HAS_RICH else None
PT_DIALOG_STYLE = PTStyle.from_dict({
    'dialog': 'bg:#000000 #d7dce2',
    'dialog frame.label': 'bg:#000000 #d7dce2 bold',
    'dialog.body': 'bg:#000000 #d7dce2',
    'dialog shadow': 'bg:#111111',
    'button': 'bg:#343a40 #d7dce2',
    'button.focused': 'bg:#7a828c #ffffff bold',
    'dialog.body button': 'bg:#343a40 #d7dce2',
    'dialog.body button.focused': 'bg:#7a828c #ffffff bold',
    'dialog.body button-arrow': 'bg:#343a40 #d7dce2',
    'dialog.body button.focused button-arrow': 'bg:#7a828c #ffffff bold',
    'radio': '#d7dce2',
    'radio-selected': '#d7dce2',
    'checkbox': '#d7dce2',
    'checkbox-selected': '#d7dce2',
    'focused': 'bg:#4d535b #ffffff',
}) if PTStyle is not None else None

def _name_recon_dir(tag, num):
    return '%s_%04d' % (tag, num)

def _supports_color():
    return sys.stdout.isatty() and os.environ.get('TERM', 'dumb') != 'dumb'

def _style(text, key):
    if not _supports_color():
        return text
    return '%s%s%s' % (COLOR[key], text, COLOR['reset'])

def _print_banner():
    if HAS_RICH:
        CONSOLE.print(Panel.fit(
            '[bold magenta]Dragonfly Reconstruction Setup[/bold magenta]\n'
            'Interactive setup for simulation or experimental data workflows',
            border_style='cyan',
        ))
        return
    print(_style('=' * 80, 'section'))
    print(_style('Dragonfly Reconstruction Setup', 'title'))
    print('Interactive setup for simulation or experimental data workflows')
    print(_style('=' * 80, 'section'))

def _print_section(title):
    if HAS_RICH:
        CONSOLE.rule('[bold cyan]%s[/bold cyan]' % title)
        return
    print()
    print(_style(title, 'section'))
    print(_style('-' * len(title), 'section'))

def _print_message(message, level=None):
    if HAS_RICH:
        if level == 'warning':
            CONSOLE.print(message, style='yellow', markup=False)
        elif level == 'error':
            CONSOLE.print(message, style='red', markup=False)
        else:
            CONSOLE.print(message, markup=False)
        return
    key = {'warning': 'warn', 'error': 'error'}.get(level)
    if key is None:
        print(message)
    else:
        print(_style(message, key))

def _print_success(label, value):
    if HAS_RICH:
        CONSOLE.print(label + ': ', end='')
        CONSOLE.print(str(value), style='green', markup=False)
        return
    print('%s: %s' % (label, _style(value, 'value')))

def _first_available_num(tag, num, prefix):
    while op.exists(op.join(prefix, _name_recon_dir(tag, num))):
        num += 1
    return num

def _legacy_create_new_recon_dir(tag='recon', num=1, prefix='./'):
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

def _parse_yes_no(raw, default=True):
    if not raw.strip():
        return default
    answer = raw.strip().lower()
    if answer in ('y', 'yes'):
        return True
    if answer in ('n', 'no'):
        return False
    raise ValueError

def _prompt_text(label, default=None, allow_empty=False):
    while True:
        if HAS_PROMPT_TOOLKIT:
            prompt_text = label
            if default is not None:
                prompt_text += ' [%s]' % default
            prompt_text += ': '
            raw = pt_prompt(prompt_text, default='' if default is None else str(default))
        else:
            prompt = label
            if default is not None:
                prompt += ' [%s]' % default
            prompt += ': '
            raw = input(_style(prompt, 'prompt'))
        if not raw.strip():
            if default is not None:
                return default
            if allow_empty:
                return ''
        elif raw.strip() or allow_empty:
            return raw.strip()
        _print_message('A value is required.', level='error')

def _prompt_path_text(label, default=None, only_directories=False):
    completer = None
    if HAS_PROMPT_TOOLKIT:
        completer = PathCompleter(only_directories=only_directories, expanduser=True)
    while True:
        if HAS_PROMPT_TOOLKIT:
            prompt_text = label
            if default is not None:
                prompt_text += ' [%s]' % default
            prompt_text += ': '
            raw = pt_prompt(prompt_text, default='' if default is None else str(default), completer=completer)
        else:
            raw = _prompt_text(label, default=default)
        if not raw.strip() and default is not None:
            return default
        if raw.strip():
            return raw.strip()
        _print_message('A value is required.', level='error')

def _prompt_yes_no(label, default=True):
    if HAS_PROMPT_TOOLKIT:
        result = yes_no_dialog(
            title='Dragonfly Setup',
            text=label,
            yes_text='Yes',
            no_text='No',
            style=PT_DIALOG_STYLE,
        ).run()
        if result is None:
            return default
        return result
    suffix = 'Y/n' if default else 'y/N'
    while True:
        raw = input(_style('%s [%s]: ' % (label, suffix), 'prompt'))
        try:
            return _parse_yes_no(raw, default=default)
        except ValueError:
            _print_message("Please respond with 'y' or 'n'.", level='error')

def _prompt_choice(label, options, default=None):
    if HAS_PROMPT_TOOLKIT:
        values = [(idx, option) for idx, option in enumerate(options, start=1)]
        result = radiolist_dialog(
            title='Dragonfly Setup',
            text=label,
            values=values,
            default=default,
            ok_text='Continue',
            cancel_text='Cancel',
            style=PT_DIALOG_STYLE,
        ).run()
        if result is None:
            raise KeyboardInterrupt
        return result
    while True:
        print(_style(label, 'prompt'))
        for idx, option in enumerate(options, start=1):
            marker = ''
            if default == idx:
                marker = ' (default)'
            print('  %d. %s%s' % (idx, option, marker))
        raw = input(_style('Choose an option: ', 'prompt')).strip()
        if not raw and default is not None:
            return default
        if raw.isdigit() and 1 <= int(raw) <= len(options):
            return int(raw)
        _print_message('Enter one of the numbered options.', level='error')

def _prompt_int(label, default=None, minimum=None):
    while True:
        raw = _prompt_text(label, default=None if default is None else str(default))
        try:
            value = int(raw)
        except ValueError:
            _print_message('Enter an integer value.', level='error')
            continue
        if minimum is not None and value < minimum:
            _print_message('Value must be at least %d.' % minimum, level='error')
            continue
        return value

def _prompt_float(label, default=None, minimum=None):
    while True:
        raw = _prompt_text(label, default=None if default is None else str(default))
        try:
            value = float(raw)
        except ValueError:
            _print_message('Enter a numeric value.', level='error')
            continue
        if minimum is not None and value < minimum:
            _print_message('Value must be at least %s.' % minimum, level='error')
            continue
        return raw

def _prompt_detsize(default='150'):
    while True:
        raw = _prompt_text('Detector size in pixels (one value or X Y)', default=default)
        parts = raw.split()
        if len(parts) not in (1, 2):
            _print_message('Enter one integer or two integers separated by spaces.', level='error')
            continue
        try:
            values = [int(part) for part in parts]
        except ValueError:
            _print_message('Detector size must contain integers only.', level='error')
            continue
        if min(values) <= 0:
            _print_message('Detector size values must be positive.', level='error')
            continue
        return ' '.join(str(value) for value in values)

def _prompt_beta_schedule(default='2.0 10'):
    jump_default, period_default = default.split()
    period = _prompt_int('Change beta_factor every how many iterations?',
                         default=int(period_default), minimum=1)
    jump = _prompt_float('Multiply beta_factor by how much at each change?',
                         default=jump_default, minimum=0)
    return '%s %d' % (jump, period)

def _prompt_recon_type(default='3d'):
    default_choice = 1 if default == '3d' else 2
    choice = _prompt_choice('Reconstruction type', ['3D', '2D'], default=default_choice)
    return '3d' if choice == 1 else '2d'

def _prompt_recon_shape_params(emc, defaults_2d=None):
    recon_type = _prompt_recon_type(default=emc.get('recon_type', '3d'))
    num_div_default = emc.get('num_div', '6')
    num_rot_default = emc.get('num_rot', RECON_DEFAULTS_2D['num_rot'])
    num_modes_default = emc.get('num_modes', RECON_DEFAULTS_2D['num_modes'])
    emc['recon_type'] = recon_type
    emc.pop('num_div', None)
    emc.pop('num_rot', None)
    emc.pop('num_modes', None)
    if recon_type == '3d':
        emc['num_div'] = str(_prompt_int('Quaternion sampling num_div', default=int(num_div_default), minimum=1))
    else:
        rot_default = num_rot_default if defaults_2d is None else defaults_2d['num_rot']
        modes_default = num_modes_default if defaults_2d is None else defaults_2d['num_modes']
        emc['num_rot'] = str(_prompt_int('Number of in-plane rotations', default=int(rot_default), minimum=1))
        emc['num_modes'] = str(_prompt_int('Number of modes', default=int(modes_default), minimum=1))

def _prompt_existing_path(label, default=None):
    while True:
        raw = _prompt_path_text(label, default=default)
        if op.exists(raw):
            return op.realpath(raw)
        _print_message('Path does not exist: %s' % raw, level='error')

def _prompt_existing_dir(label, default=None):
    while True:
        raw = _prompt_path_text(label, default=default, only_directories=True)
        if op.isdir(raw):
            return op.realpath(raw)
        _print_message('Directory does not exist: %s' % raw, level='error')

def _parse_index_selection(raw, max_index):
    selected = []
    seen = set()
    for item in raw.split(','):
        token = item.strip()
        if not token:
            continue
        if '-' in token:
            bounds = token.split('-', 1)
            if len(bounds) != 2 or not bounds[0].isdigit() or not bounds[1].isdigit():
                raise ValueError
            start = int(bounds[0])
            stop = int(bounds[1])
            if start > stop:
                raise ValueError
            values = range(start, stop + 1)
        else:
            if not token.isdigit():
                raise ValueError
            values = [int(token)]
        for value in values:
            if value < 1 or value > max_index:
                raise ValueError
            if value not in seen:
                selected.append(value)
                seen.add(value)
    if not selected:
        raise ValueError
    return selected

def _pick_files_from_directory():
    while True:
        folder = op.realpath(_prompt_existing_dir('Folder containing photon files', default='.'))
        entries = [name for name in sorted(os.listdir(folder))
                   if op.isfile(op.join(folder, name))]
        if not entries:
            _print_message('No files found in %s' % folder, level='warning')
            continue
        if HAS_PROMPT_TOOLKIT:
            values = [(name, name) for name in entries]
            selection = checkboxlist_dialog(
                title='Photon File Selection',
                text='Select files to include from %s' % folder,
                values=values,
                ok_text='Use Selection',
                cancel_text='Back',
                style=PT_DIALOG_STYLE,
            ).run()
            if selection:
                return folder, list(selection)
            _print_message('Select at least one file.', level='warning')
            continue
        _print_section('File Picker')
        print('Select files using syntax like 1,3,8-12')
        for idx, name in enumerate(entries, start=1):
            print('  %3d. %s' % (idx, name))
        while True:
            raw = input(_style('Files to include: ', 'prompt')).strip()
            try:
                selection = _parse_index_selection(raw, len(entries))
            except ValueError:
                _print_message('Invalid selection. Use syntax like 1,3,8-12.', level='error')
                continue
            return folder, [entries[idx-1] for idx in selection]

def _propose_recon_location(initial_tag, initial_num, initial_prefix):
    tag = initial_tag
    prefix = initial_prefix
    run_num = _first_available_num(tag, initial_num, prefix)
    while True:
        recon_name = _name_recon_dir(tag, run_num)
        recon_dir = op.join(prefix, recon_name)
        _print_section('Reconstruction Directory')
        if HAS_RICH:
            table = Table(show_header=False, box=None)
            table.add_row('Tag', tag)
            table.add_row('Run number', str(run_num))
            table.add_row('Parent path', prefix)
            table.add_row('Proposed directory', '[green]%s[/green]' % recon_dir)
            CONSOLE.print(table)
        if not HAS_PROMPT_TOOLKIT:
            print('Proposed directory: %s' % _style(recon_dir, 'value'))
        choice = _prompt_choice(
            'How would you like to proceed?\n\nProposed directory:\n%s' % recon_dir,
            ['Accept this directory', 'Change run number', 'Change tag', 'Change parent path'],
            default=1,
        )
        if choice == 1:
            if op.exists(recon_dir):
                _print_message('That directory already exists. Pick another run number.', level='warning')
                continue
            return tag, run_num, prefix
        if choice == 2:
            run_num = _prompt_int('Run number', default=run_num, minimum=1)
            continue
        if choice == 3:
            tag = _prompt_text('Reconstruction tag', default=tag)
            run_num = _first_available_num(tag, 1, prefix)
            continue
        prefix = op.realpath(_prompt_existing_dir('Parent directory', default=prefix))
        run_num = _first_available_num(tag, 1, prefix)

def _prompt_config_style():
    choice = _prompt_choice('Config file style', ['Keep helpful comments', 'Write a clean config'], default=1)
    return choice == 1

def _prompt_simulation_config(recon_dir, parent_dir):
    config = {section: values.copy() for section, values in SIM_DEFAULTS.items()}
    install_local_pdb = op.join(parent_dir, 'aux', '4BED.pdb')
    recon_local_pdb = op.join(recon_dir, 'aux', '4BED.pdb')

    _print_section('Simulation Model')
    model_choice = _prompt_choice('Choose the structure source', ['Fetch by PDB code', 'Use a local PDB file'], default=2)
    if model_choice == 1:
        config['make_densities']['pdb_code'] = _prompt_text('PDB code', default=config['make_densities']['pdb_code'])
    else:
        pdb_file = _prompt_existing_path('Path to local PDB file', default=recon_local_pdb)
        config['make_densities'].pop('pdb_code', None)
        if op.realpath(pdb_file) == op.realpath(install_local_pdb):
            config['make_densities']['in_pdb_file'] = 'aux/4BED.pdb'
        else:
            config['make_densities']['in_pdb_file'] = pdb_file

    _print_section('Detector Geometry')
    params = config['parameters']
    params['detd'] = _prompt_float('Detector distance in mm', default=params['detd'], minimum=0)
    params['lambda'] = _prompt_float('Photon wavelength in Angstrom', default=params['lambda'], minimum=0)
    params['detsize'] = _prompt_detsize(default=params['detsize'])
    params['pixsize'] = _prompt_float('Pixel size in mm', default=params['pixsize'], minimum=0)
    params['stoprad'] = _prompt_float('Beamstop radius in pixels', default=params['stoprad'], minimum=0)
    pol_choice = _prompt_choice('Polarization correction', ['x', 'y', 'none'], default=1)
    params['polarization'] = ['x', 'y', 'none'][pol_choice-1]

    _print_section('Simulated Data')
    make_data = config['make_data']
    make_data['num_data'] = str(_prompt_int('Number of diffraction patterns', default=int(make_data['num_data']), minimum=1))
    make_data['fluence'] = _prompt_float('Incident fluence in photons/um^2', default=make_data['fluence'], minimum=0)

    _print_section('EMC Settings')
    emc = config['emc']
    _prompt_recon_shape_params(emc)
    emc['need_scaling'] = '1' if _prompt_yes_no('Enable fluence scaling', default=True) else '0'
    _print_message('beta_start[d] is computed per frame. The iteration factor is')
    _print_message('beta_factor * beta_schedule[0]**((i-1)//beta_schedule[1]).')
    emc['beta_factor'] = _prompt_float('Initial beta_factor', default=emc['beta_factor'], minimum=0)
    emc['beta_schedule'] = _prompt_beta_schedule(default=emc['beta_schedule'])
    return config

def _prompt_experimental_photons(recon_dir):
    _print_section('Photon Inputs')
    choice = _prompt_choice(
        'How should the photon inputs be configured?',
        ['Use a single existing file', 'Use an existing list file', 'Create a new list file from a folder'],
        default=3,
    )
    if choice == 1:
        return {'in_photons_file': _prompt_existing_path('Photon file path')}
    if choice == 2:
        return {'in_photons_list': _prompt_existing_path('Photon list file path')}

    folder, files = _pick_files_from_directory()
    list_name = _prompt_text('Photon list filename', default='photons.lst')
    list_path = op.join(recon_dir, list_name)
    if op.exists(list_path) and not _prompt_yes_no('Overwrite %s?' % list_name, default=False):
        return _prompt_experimental_photons(recon_dir)
    with open(list_path, 'w') as fptr:
        for name in files:
            fptr.write(op.join(folder, name) + '\n')
    _print_success('Wrote photon list', list_path)
    return {'in_photons_list': list_name}

def _prompt_experimental_config(recon_dir):
    config = {section: values.copy() for section, values in EXP_DEFAULTS.items()}
    config['emc'].update(_prompt_experimental_photons(recon_dir))

    _print_section('Experimental EMC Settings')
    emc = config['emc']
    emc['in_detector_file'] = _prompt_existing_path('Detector file path')
    emc['output_folder'] = _prompt_text('Output folder', default=emc['output_folder'])
    emc['log_file'] = _prompt_text('Log file', default=emc['log_file'])
    _prompt_recon_shape_params(emc, defaults_2d=RECON_DEFAULTS_2D)
    emc['need_scaling'] = '1' if _prompt_yes_no('Enable fluence scaling', default=True) else '0'
    _print_message('beta_start[d] is computed per frame. The iteration factor is')
    _print_message('beta_factor * beta_schedule[0]**((i-1)//beta_schedule[1]).')
    emc['beta_factor'] = _prompt_float('Initial beta_factor', default=emc['beta_factor'], minimum=0)
    emc['beta_schedule'] = _prompt_beta_schedule(default=emc['beta_schedule'])
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
    config['emc'] = ordered_emc
    return config

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
    new_recon_dir = _legacy_create_new_recon_dir(tag=args.recon_tag, num=args.run_tag,
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
    _print_banner()
    tag, run_num, prefix = _propose_recon_location(args.recon_tag, args.run_tag, args.recon_prefix)
    keep_comments = _prompt_config_style()
    workflow = _prompt_choice('Choose the workflow', ['Simulation', 'Experimental'], default=1)

    recon_dir = _legacy_create_new_recon_dir(tag=tag, num=run_num, prefix=prefix)
    _setup_aux_dir(recon_dir, parent_dir, copy_aux=args.copy_aux)

    if workflow == 1:
        config = _prompt_simulation_config(recon_dir, parent_dir)
        workflow_name = 'simulation'
    else:
        config = _prompt_experimental_config(recon_dir)
        workflow_name = 'experimental'

    _write_generated_config(recon_dir, _render_config(config, keep_comments, workflow_name))

    _print_section('Setup Complete')
    _print_success('Created new directory', recon_dir)
    _print_success('Config file', op.join(recon_dir, 'config.ini'))
    _print_success('Next step', 'cd %s' % recon_dir)

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
