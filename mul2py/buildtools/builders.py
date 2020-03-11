import hyperspy.api as hs
import numpy as np
import h5py as hdf
from pathlib import Path


def resolve_reference(hdf_file, ref):
    """Resolve a HDF reference"""
    if isinstance(ref, hdf.h5r.Reference):
        return hdf_file[ref]
    else:
        return ref


def unpack_grp(hdf_group, hdf_file=None):
    """Return the members and sub-members of a hdf group as a dictionary

    Iteratively unpacks all the members and sub-members of a hdf group and returns it as a dictionary.
    If the HDF file object is also given, the refernces will be resolved before being added to the dictionary.

    Parameters
    ----------
    hdf_group : A h5py._hl.group.Group object
        The group to unpack
    hdf_file : A h5py file object, optional.
        If given, also attempt to resolve any references found in the hdf group.

    Returns
    -------
    dictionary : A dict of the members in the hdf5 group

    """
    if isinstance(hdf_group, hdf._hl.group.Group):
        dictionary = {}
        for key in hdf_group.keys():
            dictionary[key] = unpack_grp(hdf_group[key])
        return dictionary
    else:
        value = np.array(hdf_group)
        if hdf_file is not None:  # Attempt to resolve eventual references
            value = np.array([resolve_reference(hdf_file, ref) for ref in value.flat]).reshape(np.shape(value))

        # Extract different data shapes (so that integers are not 2D arrays etc
        if value.shape == (1, 1):
            value = value.ravel()[0]
        elif value.shape[0] > 1 and value.shape[1] == 1:
            value = value.ravel()
        else:
            value = value
        return value


def set_original_metadata(signal, dictionary, filepath='', simulation_type=''):
    try:
        signal.original_metadata.add_dictionary({'SimulationParameters': dictionary})
        signal.original_metadata.General.original_filename = filepath
        signal.original_metadata.Signal.simulation_type = simulation_type
    except KeyError as e:
        print(e)


def set_axes(signal, ax, values, name=None, units='Å', set_offset=True):
    if name is not None:
        signal.axes_manager[ax].name = name
    try:
        signal.axes_manager[ax].units = units
        if np.size(values) > 1:
            if set_offset:
                signal.axes_manager[ax].offset = np.min(values)
            signal.axes_manager[ax].scale = (np.max(values) - np.min(values)) / np.size(values)
        else:
            signal.axes_manager[ax].offset = values
            signal.axes_manager[ax].scale = 1
            signal.axes_manager[ax].units = 'px'
    except ValueError as e:
        print(e)


def build_ewrs(filepath):
    filepath = Path(filepath)
    """Build a HyperSpy stack from a stack of exit waves (real-space)"""
    with hdf.File(filepath, 'r') as hdf_file:
        signal = hs.signals.Signal2D(hdf_file['results']['images'])
        xs = np.array(hdf_file['results']['xs']).ravel()
        ys = np.array(hdf_file['results']['ys']).ravel()
        dx = np.array(hdf_file['results']['dx']).ravel()[0]
        dy = np.array(hdf_file['results']['dy']).ravel()[0]
        input_parameters = unpack_grp(hdf_file['results']['input'])
        hdf_file.close()
    signal = signal.transpose(navigation_axes=[3, 4, 0], signal_axes=[1, 2], optimize=True)

    # Set original metadata
    set_original_metadata(signal, input_parameters, filepath=filepath, simulation_type='EWRS')

    # Set axes properties
    set_axes(signal, 0, xs, name='x')
    set_axes(signal, 1, 'y', ys)
    set_axes(signal, 2, 'z', signal.original_metadata.SimulationParameters.thick)
    set_axes(signal, 3, 'dx', dx, set_offset=False)
    set_axes(signal, 4, 'dy', dy, set_offset=False)

    return signal


def build_cbed(filepath):
    filepath = Path(filepath)
    """Build a HyperSpy stack from SCBED data"""
    with hdf.File(filepath, 'r') as hdf_file:
        signal = hs.signals.Signal2D(hdf_file['results']['images'])
        x = np.array(hdf_file['results']['x']).ravel()[0]
        y = np.array(hdf_file['results']['y']).ravel()[0]
        dx = np.array(hdf_file['results']['dx']).ravel()[0]
        dy = np.array(hdf_file['results']['dy']).ravel()[0]
        input_parameters = unpack_grp(hdf_file['results']['input'])
        input_parameters['x'] = x
        input_parameters['y'] = y
        hdf_file.close()

    signal = signal.transpose(navigation_axes=[1], signal_axes=[0, 2], optimize=True)

    # Set original metadata
    set_original_metadata(signal, input_parameters, filepath=filepath, simulation_type='CBED')

    # Set axes properties
    set_axes(signal, 0, 'z', signal.original_metadata.SimulationParameters.thick)
    set_axes(signal, 1, 'x', dx)
    set_axes(signal, 2, 'y', dy)

    return signal


def build_scbed(filepath):
    filepath = Path(filepath)
    """Build a HyperSpy stack from SCBED data"""
    with hdf.File(filepath, 'r') as hdf_file:
        signal = hs.signals.Signal2D(hdf_file['results']['images'])
        xs = np.array(hdf_file['results']['xs']).ravel()
        ys = np.array(hdf_file['results']['ys']).ravel()
        dx = np.array(hdf_file['results']['dx']).ravel()[0]
        dy = np.array(hdf_file['results']['dy']).ravel()[0]
        input_parameters = unpack_grp(hdf_file['results']['input'])
        hdf_file.close()
    signal = signal.transpose(navigation_axes=[3, 4, 0], signal_axes=[1, 2], optimize=True)

    # Set original metadata

    # Set axes properties
    try:
        set_axes(signal, 0, 'x', xs)
        set_axes(signal, 1, 'y', ys)
        set_axes(signal, 2, 'z', signal.original_metadata.SimulationParameters.thick)
        set_axes(signal, 3, 'dx', dx, units='Å^-1')
        set_axes(signal, 4, 'dy', dy, units='Å^-1')
    except ValueError as e:
        print(e)

    return signal


def build_hrtem(filepath):
    filepath = Path(filepath)
    """Build a HyperSpy stack from HRTEM data"""
    with hdf.File(filepath, 'r') as hdf_file:
        signal = hs.signals.Signal2D(hdf_file['results']['images'])
        dx = np.array(hdf_file['results']['dx']).ravel()[0]
        dy = np.array(hdf_file['results']['dy']).ravel()[0]
        input_parameters = unpack_grp(hdf_file['results']['input'])
        hdf_file.close()
    signal = signal.transpose(navigation_axes=[2], signal_axes=[0, 1], optimize=True)

    # Set original metadata
    set_original_metadata(signal, input_parameters, filepath=filepath, simulation_type='HRTEM')

    # Set axes properties
    set_axes(signal, 0, 'x', dx)
    set_axes(signal, 1, 'y', dy)
    set_axes(signal, 2, 'z', signal.original_metadata.SimulationParameters.thick)

    return signal


def build_stem(filepath):
    return NotImplementedError('Conversion of STEM data is not implemented yet')


def build_ped(filepath):
    return NotImplementedError('Conversion of PED data is not implemented yet')


def build_sped(filepath):
    return NotImplementedError('Conversion of SPED data is not implemented yet')
