import numpy as np
import ase
from pathlib import Path
from scipy.constants import pi
import scipy.io as sio


def save_multem_model(filename, model, *args, dz=None, rms3d={}, B={}, **kwargs):
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

    # Calculate rms3d values from Debye-Waller factors if no rms3d values are provided
    if len(rms3d) == 0 and len(B) > 0:
        for atomic_number in B:
            rms3d[atomic_number] = np.sqrt(B[atomic_number] / (8 * pi ** 2))

    # Check that rms3D values are provided for all atom species in the model
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

    # Set sample dimensions (used in file header)
    order = kwargs.get('order', [0, 1, 2])
    lx = model.get_cell_lengths_and_angles()[order[0]]
    ly = model.get_cell_lengths_and_angles()[order[1]]
    lz = model.get_cell_lengths_and_angles()[order[2]]

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
    if lz % dz > 1E-10:
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
