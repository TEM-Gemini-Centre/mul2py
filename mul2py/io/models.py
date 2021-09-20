import numpy as np
import ase
from pathlib import Path
from scipy.constants import pi
import scipy.io as sio
from tabulate import tabulate


def save_multem_model(filename, model, *args, dz=None, rms3d={}, B={}, lx=None, ly=None, lz=None, **kwargs):
    '''Save a crystal supercell in a MULTEM-readable format

    This function loads a crystal file (in an ASE-readable format) and stores it as a MULTEM-readable format. MULTEM requires that each atom in the crystal have an RMS value (Related to the Debye-Waller factor), an occupancy, a label, and a charge. This function will append these values to each atom, where the occupancy, the label, and the chaarge is set to 1, 0, and 0, respectively. This may be ammended in the future. Regarding the RMS3D values, these are related to the temperature and should be looked up in tables. They are related to the Debye-Waller factors, and either may be proveded as keyword parameters. In addition to adding these values to the atom positions, this function also builds the header of the file to give it sensible values.

    To use the output ".mat" files from this function, please see the suggestion in the Notes section.

    Parameters
    ----------
    filename : str
        Path to store the MULTEM file
    model : ase.Atoms
        The model to be converted
    rms3d : dict, optional.
        Dictionary of atomic number and RMS value-pairs.
    B : dict, optional.
        Dictionary of atomic number and Debye-Waller factor value-pairs. Only used if no `rms3d` is given.
    lx, ly, lz: float, optional.
        The simulation box dimensions in Å. Defaults to the supercell dimensions.
    a, b, c : float, optional.
        The crystal unit cell vector lengths. Defaults to the supercell dimensions.
    na, nb, nc : float, optional.
        The number of unit cells in a, b, and c directions. Defaults to relative box lengths
    slices : float, optional.
        The number of slices. Should be larger or equal to the number of unit cells in the Z direction. Defaults to the supercell Z-dimension divided by `c`.

    Other Parameters
    ----------------
    args : optional positional arguments.
        specify 'gui' to output file for use with the MULTEM GUI.
    kwargs: optional keyword arguments.
        Specify order = [0, 1, 2] to adjust order of model size dimension if model cell dimensions are swapped in ASE (can happen when model is rotated).
        Specify 'a', 'b', or 'c' to set unit cell dimension of metadata
        Specify 'na', 'nb', or 'nc' to set number of unit cells in a, b, and c in metadata

    Notes
    -----
    Requires the ASE Atomic Simulation Environment Python package.

    To load the output of this function into MULTEM, use the following lines in the matlab script file:
    `input_multislice = multe_default_values(); model_name = dir('output_name.mat'); load(model_name.name);`
    `input_multislice.spec_atoms = spec_atoms;`
    `input_multislice.spec_lx = spec_lx;`
    `input_multislice.spec_ly = spec_ly;`
    `input_multislice.spec_lz = spec_lz;`
    `input_multislice.spec_dz = spec_dz;`

    Courtesy of Thomas Årholt. Changes by Emil Christiansen
    '''

    if isinstance(model, str):
        model = ase.io.read(Path(model))

    # This function only works with full occupancy as of yet, but should be simple to implement if required.
    occ = [1]
    # This function only works with a single label as of yet, but should be simple to implement if required.
    label = [0]
    # This function only works with charge-less atoms as of yet, but should be simple to implement if required.
    charge = [0]

    print(B)
    # Calculate rms3d values from Debye-Waller factors if no rms3d values are provided
    if len(rms3d) == 0 and len(B) > 0:
        for atomic_number in B:
            print(B[atomic_number])
            rms3d[atomic_number] = np.sqrt(B[atomic_number] / (8 * pi ** 2))

    # Check that rms3D values are provided for all atom species in the model
    print(rms3d)
    for atomic_number in set(model.get_atomic_numbers()):
        if not atomic_number in rms3d: raise ValueError(
            'No rms3d value for atoms with Z={Z:.0f}'.format(Z=atomic_number))

    # The data structure to be written to the MULTEMreadable file:
    try:
        data = [[Z] + list(xyz) + [rms3d[Z]] + occ + label + charge for Z, xyz in zip(model.numbers, model.positions)]
    except KeyError as e:
        print(
            'Could not build data structure for MULTEM file with atoms {model.numbers!r} at {model.positions!r}. Missing rms3D'.format(
                model=model))
        raise e

    order = kwargs.get('order', [0, 1, 2])
    if lx is None:
        lx = model.get_cell_lengths_and_angles()[order[0]]
    if ly is None:
        ly = model.get_cell_lengths_and_angles()[order[1]]
    if lz is None:
        lz = model.get_cell_lengths_and_angles()[order[2]]
    # Set sample dimensions (used in file header)
    print('Model size\n\tX: {:.2f} Å\n\tY: {:.2f} Å\n\tZ: {:.2f} Å'.format(lx, ly, lz))

    # Check sample dimensions
    x_max, y_max, z_max = np.max([atom.position for atom in model], axis=0)
    x_min, y_min, z_min = np.min([atom.position for atom in model], axis=0)
    message = []
    ok = False
    if abs(x_max-x_min) >= lx:
        message.append('Span {:.2f} Å between atoms in x-dimension exceeds model dimension {:.2f} Å.'.format(x_max-x_min, lx))
    elif abs(y_max-y_min) >= ly:
        message.append(
            'Span {:.2f} Å between atoms in y-dimension exceeds model dimension {:.2f} Å.'.format(y_max - y_min, ly))
    elif abs(z_max-z_min) >= lz:
        message.append(
            'Span {:.2f} Å between atoms in z-dimension exceeds model dimension {:.2f} Å.'.format(z_max - z_min, lz))
    else:
        message.append('Atom positions lies within simulation dimensions: OK')
        ok = True
    message = '\n'.join(message)
    if not ok:
        raise ValueError(message)
    else:
        print(message)

    # Get crystal unit cell vectors. Defaults to the sample dimensions. I don't think these are used directly in MULTEM anywhere.
    a = kwargs.get('a', lx)
    b = kwargs.get('b', ly)
    c = kwargs.get('c', lz)
    na = kwargs.get('na', lx / a)
    nb = kwargs.get('nb', ly / b)
    nc = kwargs.get('nc', lz / c)

    # Get the slice thickness
    if dz is None:
        slices = kwargs.get('slices', nc)
        dz = lz / slices
    if (abs(int(lz/dz) - lz/dz) > 1E-10):
        print('Warning: The slice thickness does not match the specimen z-size: lz/dz = {}!'.format(lz / dz))

    filename = Path(filename)
    # Get the output type
    if 'gui' not in args:  # Build a dictionary to save to the ".mat" file:
        matlabdict = {
            "spec_atoms": data,
            "spec_lx": lx,
            "spec_ly": ly,
            "spec_lz": lz,
            "spec_dz": dz,  # Preferably some fraction of the unit cell, or just a small distance
            "a": a,
            "b": b,
            "c": c,
            "na": na,
            "nb": nb,
            "nc": nc
        }
        sio.savemat(filename.with_suffix('.mat'), matlabdict)

    else:
        header = [[lx] + [ly] + [dz] + 5 * [0]]
        total = np.array(header + data)
        np.savetxt(filename.with_suffix('.txt'), total, fmt='%.8f', newline="\n", delimiter=" ")

def mat2cif(filename, **kwargs):
    '''Convert a MUTLEM model ".mat" file to ".cif" format.

    Parameters
    ----------
    filename : str
        Complete path, including extension, to the file to be converted. The converted file will be stored with the same name at the same location, but with the .cif extension

    Other parameters
    ----------------
    kwargs : Keyword arguments.
        Use this to specify the number of atom layers in the c-direction to be kept in the converted file by setting `nlayers=<int>`.

    Returns
    -------
    model : ase.Atoms object
        The model that was converted to .cif, after eventual termination along the c-direction.

    '''

    # Get the number of layers to keep in the model. Set to None by default, where all atoms are kept.
    nlayers = kwargs.get('nlayers', None)

    # Convert the filename to a Path object for easier distribution accross platforms
    filename = Path(filename)

    # Load the ".mat" file
    mat_file = sio.loadmat(filename)

    # Extract the positions and the atomic numbers of the atoms
    positions = mat_file['spec_atoms'][:, 1:4]
    numbers = mat_file['spec_atoms'][:, 0]

    # Extract the cell-dimensions of the model
    lx = mat_file['spec_lx'][0][0]
    ly = mat_file['spec_ly'][0][0]
    lz = mat_file['spec_lz'][0][0]

    # Create an ase.Atoms object from the model
    model = ase.Atoms(positions=positions, numbers=numbers, cell=np.array([lx, ly, lz]))
    if nlayers is not None:
        model = ase.build.cut(model, nlayers=nlayers)

    # Get the unique element symbols. Useful for later.
    unique_symbols = []
    for atom in model:
        if atom.symbol not in unique_symbols: unique_symbols.append(atom.symbol)

    # Open/create the new .cif file and start writing
    with open(filename.with_suffix('.cif'), 'w') as ciffile:
        # Write some info used by the .cif format first. Figured out by inspecting working .cif files
        ciffile.write('data_{filename.stem}\n'.format(filename=filename))
        ciffile.write('_cell_length_a\t{a:.3}\n'.format(a=lx))
        ciffile.write('_cell_length_b\t{b:.3}\n'.format(b=ly))
        ciffile.write('_cell_length_c\t{c:.3}\n'.format(c=lz))
        ciffile.write('_cell_angle_alpha\t90\n')
        ciffile.write('_cell_angle_beta\t90\n')
        ciffile.write('_cell_angle_gamma\t90\n')
        ciffile.write('\n')

        ciffile.write('_symmetry_space_group_name_H-M\t"P 1"\n')
        ciffile.write('_symmetry_int_tables_number\t1')
        ciffile.write('\n')

        ciffile.write('loop_\n')
        ciffile.write('_symmetry_equiv_pos_as_xyz\n')
        ciffile.write('\'x, y, z\'\n')
        ciffile.write('\n')

        ciffile.write('loop_\n')
        ciffile.write('_atom_site_label\n')
        ciffile.write('_atom_site_occupancy\n')
        ciffile.write('_atom_site_fract_x\n')
        ciffile.write('_atom_site_fract_y\n')
        ciffile.write('_atom_site_fract_z\n')
        ciffile.write('_atom_site_thermal_displace_type\n')
        ciffile.write('_atom_site_B_iso_or_equiv\n')
        ciffile.write('_atom_site_type_symbol\n')

        # loop over each atom type and add them in the file
        for symbol in unique_symbols:
            counter = 1  # used for counting the atoms of the current species that have been added to the file, used for specifying site labels
            for atom in model:  # loop over all the atoms
                if atom.symbol == symbol: ciffile.write(
                    '{symbol}{counter}\t{occ:.3f}\t{fract_x:.3f}\t{fract_y:.3f}\t{fract_z:.3f}\tBiso\t1.000\t{symbol}\n'.format(
                        symbol=symbol, counter=counter, occ=1, fract_x=atom.position[0] / lx,
                        fract_y=atom.position[1] / ly,
                        fract_z=atom.position[2] / lz))  # If the atom is of the "current" species, add it to the file
        ciffile.write('\n')  # Finish by writing an empty line to the file. Not sure if necessary...
        ciffile.close()
    return model

def verify_model(model, nx=1024, ny=1024, bwl=False, acceleration_voltage=200):
    """
    Verify the model for multisice simulations

    Check the model dimensions and outputs the maximum scattering angles and resolution a multislice image simulation of the model will have.

    Parameters
    ----------
    model: ase.Atoms
        The model to verify
    nx, ny: int, optional
        The potential sampling you will use. Defaults to 1024
    bwl: bool, optional.
        Will you bandwidth limit the simulation? Defaults to False
    acceleration_voltage: float, optional.
        The acceleration voltage of the simulation in kV. Defaults 200.

    Returns
    -------
    info: str
        The multislice info.
    """

    if bwl:
        bandwidth = 2/3
    else:
         bandwidth = 1

    wavelength = energy2wavelength(acceleration_voltage)

    Lx, Ly, Lz = np.sum(model.cell.array, axis=0)

    dx, dy = Lx/nx, Ly/ny
    dkx, dky = 1/Lx, 1/Ly
    Kx, Ky = dkx * nx * bandwidth / 2, dky * ny * bandwidth / 2
    tx, ty = Kx*wavelength, Ky*wavelength

    info = tabulate([
        ['Acceleration voltage', 'kV', '{:.0f}'.format(acceleration_voltage), '', ''],
        ['Wavelength', 'Å', '{:.4f}'.format(wavelength), '', ''],
        ['Bandwidth', '', '{:.2f}'.format(bandwidth), '', ''],
        ['Model size', 'Å', '', '{:.4f}'.format(Lx), '{:.4f}'.format(Ly)],
        ['Potential sampling', 'px', '', '{:.0f}'.format(nx), '{:.0f}'.format(ny)],
        ['Real-space resolution', 'Å', '', '{:.4f}'.format(dx), '{:.4f}'.format(dy)],
        ['Reciprocal resolution', '1/Å', '', '{:.4f}'.format(dkx), '{:.4f}'.format(dky)],
        ['Maximum scattering vector', '1/Å', '', '{:.4f}'.format(Kx), '{:.4f}'.format(Ky)],
        ['Maximum scattering angle', 'mrad', '', '{:.4f}'.format(tx), '{:.4f}'.format(ty)]
    ],
        headers = ['', 'Unit', '', 'X', 'Y']
    )

    return info

def make_model(unit_cell, zone_axis, upwards_direction, nx=1, ny=1, nz=1):
    """
    Make a model from a unit cell.

    Creates a model with the given zone axis and upwards direction.

    Parameters
    ----------
    unit_cell: ase.atoms or str
        An ASE atoms object or a path to a unit cell file.
    zone_axis: list, tuple, array
        The zone axis of the model given as a 1D array-like object of length 3
    upwards_direction: list, tuple, array
        The upwards direction of the model given as a 1D array-like object of length 3. Should be perpendicular to the zone axis.
    nx, ny, nz: int, optional
        The number of lattice planes (given by the zone axis, upwards direction, and their cross-product) in x, y, and z directions. Defaults to 1.

    Returns
    -------
    model: ase.Atoms
        The model
    size: 3-tuple
        The size of the model in x, y, and z measured in Å.
    """

    nx=int(nx)
    ny=int(ny)
    nz=int(nz)

    zone_axis=np.array(zone_axis)
    upwards_direction=np.array(upwards_direction)

    #Check dimensions of vectors
    if not np.shape(zone_axis) == (3,):
        raise ValueError('Zone axis {P} must have shape (3,), got {s}'.format(P=zone_axis, s=np.shape(zone_axis)))

    if not np.shape(upwards_direction) == (3,):
        raise ValueError('Upward direction {P} must have shape (3,), got {s}'.format(P=zone_axis, s=np.shape(upwards_direction)))

    #Check perpendicularity
    if not abs(np.dot(zone_axis, upwards_direction)) <1E-10:
        raise ValueError('Zone axis and upwards direction are not perpendicular (dot product: {d:.12e}'.format(d=np.dot(zone_axis, upwards_direction)))

    in_plane_direction = np.cross(upwards_direction, zone_axis)

    if not isinstance(unit_cell, ase.Atoms):
        unit_cell = ase.io.read(unit_cell)

    print('Making model from {unit_cell} with:'.format(unit_cell=unit_cell))
    table = tabulate([['Zone axis', zone_axis[0], zone_axis[1], zone_axis[2]],
                      ['Upwards direction', upwards_direction[0], upwards_direction[1], upwards_direction[2]],
                      ['In-plane direction', in_plane_direction[0], in_plane_direction[1], in_plane_direction[2]],
                      ['Size', nx, ny, nz]], headers=['', 'X', 'Y', 'Z'])
    print(table)

    model = cut(unit_cell, a=in_plane_direction*nx, b=zone_axis*nz, c=upwards_direction*ny) #Cut unit cell into given size and shape

    model.rotate(zone_axis, 'z', rotate_cell=True) #Rotate the zone axis into the Z direction
    model.rotate(upwards_direction, 'y', rotate_cell=True) #Rotate the upwards-direction into the Y direction
    model.center() #Center the atoms

    lx, ly, lz = np.sum(model.cell.array, axis=0) #Get the model dimensions. For an _incorrect_ model where the model edges are not aligned with X, Y, or Z, some of these might become zero.

    #Print some info
    print('Made model:')
    x_max, y_max, z_max = np.max([atom.position for atom in model], axis=0) #max atomic position#Print some info
    print(tabulate([['Model size [Å]', lx, ly, lz], ['Max atom pos [Å]', x_max, y_max, z_max]], headers = ['', 'X', 'Y', 'Z']))

    return model, (lx, ly, lz)

def energy2wavelength(tension, e=1.602E-19, m0=9.109E-31, h=6.626E-34, c=2.998E8, **kwargs):
    """Convert acceleration voltage to electron wavelength

    Parameters
    ---------
    tension: float
        Acceleration voltage in kV
    e: float, optional
        Electron charge, defaults to 1.602E-19 C
    m0: float, optional
        Rest mass of electron, defaults to 9.109E-31 kg
    h: float, optional
        Plancks constant, defaults to 6.626E-34 N m s
    c: float, optional
        Speed of light in vacuum, defaults to 2.998E8 m/s

    Other Parameters
    ----------------
    **kwargs: Optional keyword arguments
        unit: str.
            The unit of the returned wavelength, should be 'm', 'nm', 'Å', or 'pm'. Default is 'Å'

    Returns
    -------
    wavelength: float
        The electron wavelength in [m], [nm], [Å] (default), or [pm], depending on the value passed to kwargs['unit'].
    """
    tension = tension * 1E3 #kV to V

    unit = str(kwargs.get('unit', 'Å'))
    conversion_factors = {
        'm': 1,
        'nm': 1E9,
        'Å': 1E10,
        'pm': 1E12
    }

    return h / (np.sqrt(2 * m0 * e * tension * (1. + e * tension / (2. * m0 * c**2)))) * conversion_factors[unit]