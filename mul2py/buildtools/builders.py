import hyperspy.api as hs
import numpy as np
from pathlib import Path
from functools import reduce
from ..io.hdf import HDFReader


class Error(Exception):
    pass


class LoadError(Error):
    pass


def set_original_metadata(signal, dictionary, filepath='', simulation_type='', elapsed_time=None):
    filepath = Path(filepath)
    try:
        signal.original_metadata.add_dictionary({'SimulationParameters': dictionary})
        signal.original_metadata.add_dictionary({
            'General': {
                'original_filename': str(filepath),
                'elapsed_time': elapsed_time,
                'title': str(filepath.stem)
            }
        })
        signal.original_metadata.add_dictionary({'Signal': {'simulation_type': simulation_type}})
    except KeyError as e:
        print(e)


def set_important_simulation_parameters(signal, simulation_type=None, metadata_dict=None):
    """Sets the signal metadata field `SimulationParameters` based on contents in the original metadata

    Goes through predefined metadata keys that are deemed important for the given simulation type and copies the relevant metadata from `signal.original_metadata`.

    Parameters
    ----------
    signal: hyperspy.api.signals.Signal2D
        The signal to set metadata values for
    simulation_type: str or None. Optional.
        The simulation type. THis will decide which parameters are set. Default is None, in which case it sets "common" important parameters such as the potential grid resolution and specimen parameters.
    metadata_dict: dict or None. Optional.
        A dictionary with key-value pairs for manual specification of metadata. Can be used to specify user-specific details, such as the GPU that was used, the simulation time, etc.


    """
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
            'iw_x',
            'iw_y'
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

    if simulation_type is None:
        simulation_type = 'all'

    if metadata_dict is not None:
        for key in metadata_dict:
            try:
                signal.metadata.add_dictionary({'SimulationParameters': {key: metadata_dict[key]}})
            except KeyError as e:
                print('Exception when setting manual simulation parameter "{key}". \n{error}'.format(key=key, error=e))
    else:
        for important_parameter in _important_parameters[simulation_type]:
            try:
                signal.metadata.add_dictionary({
                    'SimulationParameters': {
                        important_parameter: signal.original_metadata['SimulationParameters'][important_parameter]
                    }
                })
            except KeyError as e:
                print(
                    'Exception when setting Simulation parameter metadata for important parameter "{key}":\n{error}'.format(
                        key=important_parameter, error=e))
                pass


def set_axes(signal, ax, values, **kwargs):
    name = kwargs.get('name', '')
    units = kwargs.get('units', 'Å')
    offset = kwargs.get('offset', None)
    scale = kwargs.get('scale', None)

    if np.size(values) > 1:
        if offset is None:
            offset = np.min(values)
        if scale is None:
            scale = (np.max(values) - np.min(values)) / np.size(values)
    else:
        if offset is None:
            offset = values
        if scale is None:
            scale = values

    signal.axes_manager[ax].name = name
    signal.axes_manager[ax].units = units
    signal.axes_manager[ax].scale = scale
    signal.axes_manager[ax].offset = offset


def load_output(output_file, image_type):
    """ Reads a MULTEM `output_multislice` file

    Reads a hdf file and tries to build a hyperspy signal by assuming that the HDF file has a group named 'output_multislice' at its root.

    Parameters
    ----------
    output_file: str,
        Path to the HDF file
    image_type: str. 'image_tot' or 'm2psi_tot'.
        The type of data.

    Returns
    -------
    signal: hyperspy.api.signals.Signal2D
        The image stack as a hyperspy signal.
    hdf_file: HDFReader,
        The contents of the hdf file - useful for getting more details/metadata.
    """

    filepath = Path(output_file)
    hdf_file = HDFReader(filepath)
    hdf_file.remove_excess_dimensions()  # Remove the excess dimensions of all attributes

    if image_type == 'image_tot':  # Load pattern of `image_tot`
        image_shape = hdf_file.output_multislice.data.image_tot.image_tot0.resolved_reference.image.image0.resolved_reference.shape
        detectors = hdf_file.output_multislice.data.image_tot.image_tot0.resolved_reference.image.depth
        stack_depth = hdf_file.output_multislice.data.image_tot.depth
        if detectors > 1:
            stack_dimensions = (stack_depth, detectors, image_shape[0], image_shape[1])
        else:
            stack_dimensions = (stack_depth, image_shape[0], image_shape[1])
        image_stack = np.zeros(stack_dimensions)
        for t in range(stack_depth):
            for detector in range(detectors):
                # Get the image for each thickness and detector. Must use reduce and getattr to access the data properly.
                image_stack[t, detector, :, :] = reduce(getattr, (
                    'output_multislice', 'data', 'image_tot', 'image_tot{t}'.format(t=t), 'resolved_reference',
                    'image', 'image{detector}'.format(detector=detector), 'resolved_reference', '_content'),
                                                        hdf_file)
    elif image_type == 'm2psi_tot':  # Load pattern of `m2psi_tot`
        image_shape = hdf_file.output_multislice.data.m2psi_tot.m2psi_tot0.resolved_reference.shape
        stack_depth = hdf_file.output_multislice.data.m2psi_tot.depth
        stack_dimensions = (stack_depth, image_shape[0], image_shape[1])
        image_stack = np.zeros(stack_dimensions)
        for t in range(stack_depth):
            # Get the image for each thickness. Must use reduce and getattr to access the data properly.
            image_stack[t, :, :] = reduce(getattr, (
                'output_multislice', 'data', 'm2psi_tot', 'm2psi_tot{t}'.format(t=t), 'resolved_reference', '_content'),
                                          hdf_file)
    else:
        return NotImplemented
    signal = hs.signals.Signal2D(image_stack)
    hdf_file.close()
    return signal, hdf_file.output_multislice


def load_input(input_file):
    """Loads a MULTEM `input_multislice` file

    Reads a hdf file and tries to extract the contents into a dictionary by assuming that the HDF file has a group named 'input_multislice' at its root.

    Parameters
    ----------
    input_file: str
        Path to HDF file

    Returns
    -------
    input_parameters: dict
        Dictionary of the HDF file.
    """
    filepath = Path(input_file)
    hdf_file = HDFReader(filepath)
    hdf_file.remove_excess_dimensions()  # Remove the excess dimensions of all attributes

    input_parameters = hdf_file.content2dict()
    hdf_file.close()
    return input_parameters


def load_results(results_file):
    """Reads a MULTEM results file (following the ".ecmat" format)

    Reads a hdf file and builds a hyperspy signal by assumint that the HDF file has a group named 'results/images' at its root.

    The .ecmat format is my (Emil Christiansen) own way of structuring output in .mat files from MULTEM (ecmat = ECmat = easymat) that makes it easier to load e.g. scanning SCBED or scanning EWRS results.

    Parameters
    ----------
    results_file: str
        Path to HDF file.

    Returns
    -------
    signal: hs.signals.Signal2D
        The image stack as a hyperspy signal
    hdf_file: HDFReader
        The contents of the hdf file - useful for getting more details/metadata later

    Notes
    -----
    The `ecmat` format should appear like this:
    root - `results`:
        results - `images`: (nx, ny, (sx), (sy), (thick)), where nx and ny are the potential size, (sx) and (sy) are the (optional) scan dimensions and (thick) is the (optional) thickness-dimensions.
            IMPORTANT! due to differences in MATLAB/python arrays, `results.images` should contain the _transposed_ array, i.e. `results.images(:, :, i, j, t) = transpose(image);`  is the proper way of assigning in MATLAB.
        results - `input`: the `input_multislice` struct used in MULTEM (used for building metadata)
        results - `dx`: the potential resolution in x (float)
        results - `dy`: the potential resolution in y (float)
        results - `xs`: the _manual_ scan x positions (used in e.g. scanning EWRS or SCBED)
        results - `ys`: the _manual_ scan y positions (used in e.g. scanning EWRS or SCBED)
        results - `thick`: the specimen thicknesses taken from the output of the multislice simulation (optional if `results.thicknesses` is given)
        results - `thicknesses`: the specimen thicknesses taken from the output of the multislice simulation at each _manual_ scan position (for validation and completeness) (optional if `results.thick` is given)
    """

    filepath = Path(results_file)
    hdf_file = HDFReader(filepath)
    hdf_file.remove_excess_dimensions()
    image_stack = reduce(getattr, ('results', 'images', '_content'), hdf_file)
    hdf_file.close()
    return hs.signals.Signal2D(image_stack), hdf_file.results


def make_signal(filepath):
    """Loads and makes a HyperSpy signal from a '.ecmat' file"""
    _simulation_types = {
        11: 'STEM',
        12: 'ISTEM',
        21: 'SCBED',
        22: 'CBEI',
        31: 'ED',
        32: 'HRTEM',
        41: 'PED',
        42: 'HCTEM',
        51: 'EWFS',
        52: 'EWRS',
        61: 'EELS',
        62: 'EFTEM',
        71: 'ProbeFS',
        72: 'ProbeRS',
        81: 'PPFS',
        82: 'PPRS',
        91: 'TFFS',
        92: 'TFRS'
    }
    filepath = Path(filepath)
    if filepath.suffix == '.ecmat':
        signal, data = load_results(filepath)
    else:
        raise ValueError('File must be a `.ecmat` file, got `{}`'.format(filepath.suffix))

    # Get the simulation type
    try:
        simulation_type = data.simulation_type().lower()
    except AttributeError:
        try:
            simulation_type = _simulation_types[int(data.input.simulation_type())].lower()
        except KeyError as e:
            print('No simulation type for simulation key "{sim_type}":\n{err}'.format(
                sim_type=data.input.simulation_type(), err=e))
    simulation_dimension = len(np.shape(signal))

    try:
        z = data.thick()
    except AttributeError:
        z = data.thicknesses()
    dx = data.dx()
    dy = data.dy()

    # Set axes properties
    if simulation_type == 'ewrs':
        xs = data.xs()
        ys = data.ys()
        if simulation_dimension == 5:
            set_axes(signal, 0, xs, name='x')
            set_axes(signal, 1, ys, name='y')
            set_axes(signal, 2, z, name='z')
            set_axes(signal, 3, dx, name='X')
            set_axes(signal, 4, dy, name='Y')
        elif simulation_dimension == 4:
            set_axes(signal, 0, xs, name='x')
            set_axes(signal, 1, ys, name='y')
            set_axes(signal, 3, dx, name='X')
            set_axes(signal, 4, dy, name='Y')
        elif simulation_dimension == 3:
            set_axes(signal, 0, z, name='z')
            set_axes(signal, 1, dx, name='X')
            set_axes(signal, 2, dy, name='Y')
        elif simulation_dimension == 2:
            set_axes(signal, 0, dx, name='X')
            set_axes(signal, 1, dy, name='Y')
        else:
            print('Simulation dimension {:.0f} is not supported for simulation type "{}": Axes are not set.'.format(
                simulation_type))
    elif (simulation_type == 'cbed' or simulation_type == 'ped'):
        xs = data.xs()
        ys = data.ys()
        if simulation_dimension == 5:
            set_axes(signal, 0, xs, name='x')
            set_axes(signal, 1, ys, name='y')
            set_axes(signal, 2, z, name='z')
            set_axes(signal, 3, dx, name='kx', units='Å^-1')
            set_axes(signal, 4, dy, name='ky', units='Å^-1')
        elif simulation_dimension == 4:
            set_axes(signal, 0, xs, name='x')
            set_axes(signal, 1, ys, name='y')
            set_axes(signal, 2, dx, name='kx', units='Å^-1')
            set_axes(signal, 3, dy, name='ky', units='Å^-1')
        elif simulation_dimension == 3:
            set_axes(signal, 0, z, name='z')
            set_axes(signal, 1, dx, name='kx', units='Å^-1')
            set_axes(signal, 2, dy, name='ky', units='Å^-1')
        elif simulation_dimension == 2:
            set_axes(signal, 1, dx, name='kx', units='Å^-1')
            set_axes(signal, 2, dy, name='ky', units='Å^-1')
        else:
            print('Simulation dimension {:.0f} is not supported for simulation type "{}": Axes are not set.'.format(
                simulation_type))
    elif simulation_type == 'stem':
        xs = np.linspace(data.input.scanning_x0(), data.input.scanning_xe(), int(data.input.scanning_ns()),
                         endpoint=data.input.scanning_periodic())
        ys = np.linspace(data.input.scanning_y0(), data.input.scanning_ye(), int(data.input.scanning_ns()),
                         endpoint=data.input.scanning_periodic())
        if simulation_dimension == 4:
            set_axes(signal, 0, z, name='z')
            set_axes(signal, 1, 1, name='detector', scale=1, offset=1, units='')
            set_axes(signal, 2, xs, name='x')
            set_axes(signal, 3, ys, name='y')
        elif simulation_dimension == 3:
            set_axes(signal, 0, z, name='z')
            set_axes(signal, 1, xs, name='x')
            set_axes(signal, 2, ys, name='y')
        elif simulation_dimension == 2:
            set_axes(signal, 0, xs, name='x')
            set_axes(signal, 1, ys, name='y')
        else:
            print('Simulation dimension {:.0f} is not supported for simulation type "{}": Axes are not set.'.format(
                simulation_type))
    elif simulation_type == 'hrtem':
        if simulation_dimension == 3:
            set_axes(signal, 0, z, name='z')
            set_axes(signal, 1, dx, name='x')
            set_axes(signal, 2, dy, name='y')
        elif simulation_dimension == 2:
            set_axes(signal, 0, dx, name='x')
            set_axes(signal, 1, dy, name='y')
        else:
            print('Simulation dimension {:.0f} is not supported for simulation type "{}": Axes are not set.'.format(
                simulation_type))
    else:
        print('Did not recognize simulation type "{}" as a valid simulation type: Axes are not set.'.format(
            simulation_type))

    input_parameters = data.input.content2dict()

    try:
        elapsed_time = data.elapsed_time()
    except AttributeError:
        elapsed_time = None

    # Set original metadata
    set_original_metadata(signal, input_parameters, filepath=filepath, simulation_type=simulation_type,
                          elapsed_time=elapsed_time)

    # Copy the general metadata
    signal.metadata.add_dictionary({'General': signal.original_metadata.General.as_dictionary()})

    # Set important parameters (common for all simulation types)
    try:
        set_important_simulation_parameters(signal)
    except AttributeError as e:
        print('Could not set general important simulation parameter:\n{err}'.format(err=e))
    # Set simulation specific parameters
    try:
        set_important_simulation_parameters(signal, simulation_type)
    except AttributeError as e:
        print(
            'Could not set important "{sim_type}" simulation parameter:\n{err}'.format(sim_type=simulation_type, err=e))
    return signal


def build_ewrs(filepath, simulation_parameter_file=None):
    """Build a HyperSpy stack from a stack of exit waves (real-space)"""
    try:
        signal, data = load_output(filepath, 'm2psi_tot')
    except AttributeError:
        signal, data = load_results(filepath)
    finally:
        xs = data.xs()
        ys = data.ys()
        try:
            z = data.thick()
        except AttributeError:
            z = data.thicknesses()
        dx = data.dx()
        dy = data.dy()

    if simulation_parameter_file is not None:
        input_parameters = load_input(simulation_parameter_file)
    else:
        input_parameters = data.input.content2dict()

    try:
        elapsed_time = data.elapsed_time()
    except AttributeError:
        elapsed_time = None

    # Set original metadata
    set_original_metadata(signal, input_parameters, filepath=filepath, simulation_type='EWRS',
                          elapsed_time=elapsed_time)

    # Copy the general metadata
    signal.metadata.add_dictionary({'General': signal.metadata.General.as_dictionary()})

    # Set important parameters (common for all simulation types)
    try:
        set_important_simulation_parameters(signal)
    except AttributeError as e:
        print('Could not set general important simulation parameter:\n{err}'.format(err=e))
    # Set simulation specific parameters
    try:
        set_important_simulation_parameters(signal, 'ewrs')
    except AttributeError as e:
        print('Could not set important EWRS simulation parameter:\n{err}'.format(err=e))

    # Set axes properties
    set_axes(signal, 0, xs, name='x')
    set_axes(signal, 1, ys, name='y')
    set_axes(signal, 2, z, name='z')
    set_axes(signal, 3, dx, name='X')
    set_axes(signal, 4, dy, name='Y')

    return signal


def build_cbed(filepath, simulation_parameter_file=None):
    """Build a HyperSpy stack from SCBED data"""
    try:
        signal, data = load_output(filepath, 'm2psi_tot')
    except AttributeError:
        signal, data = load_results(filepath)
    finally:
        try:
            z = data.thick()
        except AttributeError:
            z = data.thicknesses()
        dx = data.dx()
        dy = data.dy()

    if simulation_parameter_file is not None:
        input_parameters = load_input(simulation_parameter_file)
    else:
        input_parameters = data.input.content2dict()

    try:
        elapsed_time = data.elapsed_time()
    except AttributeError:
        elapsed_time = None

    # Set original metadata
    set_original_metadata(signal, input_parameters, filepath=filepath, simulation_type='SCBED',
                          elapsed_time=elapsed_time)

    # Copy the general metadata
    signal.metadata.add_dictionary({'General': signal.metadata.General.as_dictionary()})

    # Set important parameters (common for all simulation types)
    try:
        set_important_simulation_parameters(signal)
    except AttributeError as e:
        print('Could not set general important simulation parameter:\n{err}'.format(err=e))
    # Set simulation specific parameters
    try:
        set_important_simulation_parameters(signal, 'cbed')
    except AttributeError as e:
        print('Could not set important SCBED simulation parameter:\n{err}'.format(err=e))

    # Set axes properties
    set_axes(signal, 0, z, name='z')
    set_axes(signal, 1, dx, name='kx', units='Å^-1')
    set_axes(signal, 2, dy, name='ky', units='Å^-1')

    return signal


def build_scbed(filepath, simulation_parameter_file=None):
    """Build a HyperSpy stack from SCBED data"""
    try:
        signal, data = load_output(filepath, 'm2psi_tot')
    except AttributeError:
        signal, data = load_results(filepath)
    finally:
        xs = data.xs()
        ys = data.ys()
        try:
            z = data.thick()
        except AttributeError:
            z = data.thicknesses()
        dx = data.dx()
        dy = data.dy()

    if simulation_parameter_file is not None:
        input_parameters = load_input(simulation_parameter_file)
    else:
        input_parameters = data.input.content2dict()

    try:
        elapsed_time = data.elapsed_time()
    except AttributeError:
        elapsed_time = None

    # Set original metadata
    set_original_metadata(signal, input_parameters, filepath=filepath, simulation_type='SCBED',
                          elapsed_time=elapsed_time)

    # Copy the general metadata
    signal.metadata.add_dictionary({'General': signal.metadata.General.as_dictionary()})

    # Set important parameters (common for all simulation types)
    try:
        set_important_simulation_parameters(signal)
    except AttributeError as e:
        print('Could not set general important simulation parameter:\n{err}'.format(err=e))
    # Set simulation specific parameters
    try:
        set_important_simulation_parameters(signal, 'scbed')
    except AttributeError as e:
        print('Could not set important SCBED simulation parameter:\n{err}'.format(err=e))

    # Set axes properties
    set_axes(signal, 0, xs, name='x')
    set_axes(signal, 1, ys, name='y')
    set_axes(signal, 2, z, name='z')
    set_axes(signal, 3, dx, name='kx', units='Å^-1')
    set_axes(signal, 4, dy, name='ky', units='Å^-1')

    return signal


def build_hrtem(filepath, simulation_parameter_file=None):
    """Build a HyperSpy stack from HRTEM data"""
    try:
        signal, data = load_output(filepath, 'm2psi_tot')
    except AttributeError:
        signal, data = load_results(filepath)
    finally:
        try:
            z = data.thick()
        except AttributeError:
            z = data.thicknesses()
        dx = data.dx()
        dy = data.dy()

    if simulation_parameter_file is not None:
        input_parameters = load_input(simulation_parameter_file)
    else:
        input_parameters = data.input.content2dict()

    try:
        elapsed_time = data.elapsed_time()
    except AttributeError:
        elapsed_time = None

    # Set original metadata
    set_original_metadata(signal, input_parameters, filepath=filepath, simulation_type='HRTEM',
                          elapsed_time=elapsed_time)

    # Copy the general metadata
    signal.metadata.add_dictionary({'General': signal.metadata.General.as_dictionary()})

    # Set important parameters (common for all simulation types)
    try:
        set_important_simulation_parameters(signal)
    except AttributeError as e:
        print('Could not set general important simulation parameter:\n{err}'.format(err=e))
    # Set simulation specific parameters
    try:
        set_important_simulation_parameters(signal, 'hrtem')
    except AttributeError as e:
        print('Could not set important HRTEM simulation parameter:\n{err}'.format(err=e))

    # Set axes properties
    set_axes(signal, 0, z, name='z')
    set_axes(signal, 1, dx, name='x')
    set_axes(signal, 2, dy, name='y')

    return signal


def build_stem(filepath, simulation_parameter_file=None):
    """Build a HyperSpy stack from STEM data"""
    try:
        signal, data = load_output(filepath, 'image_tot')
    except AttributeError:
        signal, data = load_results(filepath)
    finally:
        try:
            z = data.thick()
        except AttributeError:
            z = data.thicknesses()

        xs = np.linspace(data.input.scanning_x0(), data.input.scanning_xe(), int(data.input.scanning_ns()),
                         endpoint=data.input.scanning_periodic())
        ys = np.linspace(data.input.scanning_y0(), data.input.scanning_ye(), int(data.input.scanning_ns()),
                         endpoint=data.input.scanning_periodic())

    if simulation_parameter_file is not None:
        input_parameters = load_input(simulation_parameter_file)
    else:
        input_parameters = data.input.content2dict()

    try:
        elapsed_time = data.elapsed_time()
    except AttributeError:
        elapsed_time = None

    # Set original metadata
    set_original_metadata(signal, input_parameters, filepath=filepath, simulation_type='STEM',
                          elapsed_time=elapsed_time)

    # Copy the general metadata
    signal.metadata.add_dictionary({'General': signal.metadata.General.as_dictionary()})

    # Set important parameters (common for all simulation types)
    try:
        set_important_simulation_parameters(signal)
    except AttributeError as e:
        print('Could not set general important simulation parameter:\n{err}'.format(err=e))
    # Set simulation specific parameters
    try:
        set_important_simulation_parameters(signal, 'stem')
    except AttributeError as e:
        print('Could not set important STEM simulation parameter:\n{err}'.format(err=e))

    # Set axes properties
    set_axes(signal, 0, z, name='z')
    set_axes(signal, 1, 1, name='detector', scale=1, offset=1, units='')
    set_axes(signal, 2, xs, name='x')
    set_axes(signal, 3, ys, name='y')

    return signal


def build_ped(filepath):
    return NotImplementedError('Conversion of PED data is not implemented yet')


def build_sped(filepath):
    return NotImplementedError('Conversion of SPED data is not implemented yet')
