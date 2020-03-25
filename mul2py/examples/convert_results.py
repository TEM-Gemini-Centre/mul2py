import mul2py as m2p
import argparse
from pathlib import Path

_important_parameters = {
    'hrtem': [
        'obj_lens_c_10',
        'obj_lens_c_30',
    ],
    'stem': [
        'cond_lens_c_10',
        'cond_lens_c_30',
        'cond_lens_outer_aper_ang',
        'detector',
        'scanning_x0',
        'scanning_y0',
        'scanning_xe',
        'scanning_ye',
        'scanning_ns',
        'scanning_periodic',
    ],
    'cbed': [
        'cond_lens_c_10',
        'cond_lens_c_30',
        'cond_lens_outer_aper_ang',
        'x',
        'y'
    ],
    'scbed': [
        'cond_lens_c_10',
        'cond_lens_c_30',
        'cond_lens_outer_aper_ang',
    ],
    'ewrs': [
        'cond_lens_c_10',
        'cond_lens_c_30',
        'cond_lens_outer_aper_ang'
    ],
    'all': [
        'E_0',
        'nx',
        'ny',
        'spec_lx',
        'spec_ly',
        'spec_lz',
        'spec_dz',
        'thick',
        'thick_type',
        'spec_atoms'
    ]

}

if __name__ == '__main__':
    # Parser arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file_path', help='Filepath to the result file')
    parser.add_argument('sim_type', help='The simulation type, e.g. "SCBED", "HRTEM", or "STEM"')

    parser.add_argument('--overwrite', nargs='?', default=True,
                        help='Overwrite existing output files? Default is True')

    # Parse arguments
    arguments = parser.parse_args()
    filepath = Path(arguments.file_path)
    sim_type = str(arguments.sim_type).lower()
    sim_overwrite = bool(arguments.overwrite)

    args = vars(arguments)
    print('\nGot {} arguments:'.format(len(args)))
    for arg in args:
        print('\t{} = {}'.format(arg, args[arg]))
    print('\n')

    try:
        important_metadata_keys = _important_parameters[sim_type]
    except KeyError:
        print('No important parameters defined for simulations of type {}'.format(sim_type))
        important_metadata_keys = []

    if sim_type == 'scbed':
        signal = m2p.buildtools.builders.build_scbed(filepath)
    elif sim_type == 'hrtem':
        signal = m2p.buildtools.builders.build_hrtem(filepath)
    elif sim_type == 'stem':
        signal = m2p.buildtools.builders.build_stem(filepath)
    elif sim_type == 'cbed':
        signal = m2p.buildtools.builders.build_cbed(filepath)
    elif sim_type == 'ped':
        signal = m2p.buildtools.builders.build_ped(filepath)
    elif sim_type == 'sped':
        signal = m2p.buildtools.builders.build_sped(filepath)
    elif sim_type == 'ewrs':
        signal = m2p.buildtools.builders.build_ewrs(filepath)
    else:
        raise ValueError('Simulation type {} not recognized'.format(sim_type))

    # Set simulation specific metadata
    for important_parameter in important_metadata_keys:
        try:
            signal.metadata.add_dictionary({
                'SimulationParameters': {
                    important_parameter: signal.original_metadata['SimulationParameters'][important_parameter]}
            })
        except KeyError:
            print('Could not add parameter {key} to signal metadata, it is missing from the original metadata\n{orig_metadata!r}'.format(key=important_parameter,
                                            orig_metadata=signal.original_metadata.as_dictionary()))
        except AttributeError:
            print('Failed setting metadata for "{}", skipping...'.format(important_parameter))

    # Set general simulation metadata
    for important_parameter in _important_parameters['all']:
        try:
            signal.metadata.add_dictionary({
                'SimulationParameters': {
                    important_parameter: signal.original_metadata.as_dictionary()['SimulationParameters'][important_parameter]
                }
            })
        except KeyError:
            print('Could not add parameter {key} to signal metadata, it is missing from the original metadata\n{orig_metadata!r}'.format(key=important_parameter,
                                            orig_metadata=signal.original_metadata.as_dictionary()))
        except AttributeError:
            print('Failed setting metadata for "{}", skipping...'.format(important_parameter))

    signal.metadata.General.title = filepath.stem
    signal.metadata.Signal.simulation_type = sim_type

    signal.save(filepath.with_suffix('.hspy'), overwrite=sim_overwrite)
