from CifFile import ReadCif
from pathlib import Path
from ase.spacegroup import Spacegroup
from tabulate import tabulate
from warnings import warn
from ase import Atom
from ase.cell import Cell
from scipy.io import savemat
import numpy as np


class CIFAtom(Atom):
    def __init__(self, cell, atom=None, fractional_position=np.array([0, 0, 0]), dwf=np.nan, occupancy=1, charge=0,
                 site_label='',
                 symbol=''):
        '''
        Create a CifAtom object

        A CIFAtom is an ASE Atom object with additional fields read from a cif file, such as Debye-Waller factors and site labels.

        Parameters
        ----------
        cell: ase.cell.Cell
            The cell in which the atom belongs.
        kwargs: Keyword arguments passed on to ase.Atom
            Can also be used to specify additional fields not supported by ASE, such as Debye-Waller factors.

        Examples
        --------
        To make an atom with a debye-waller factor different form NaN (default), pass on the debye-waller factor using the keyword "debye_waller_factor":
        ```
        my_cif_atom = CIFAtom(symbol='Al', position = [0, 0, 0], debye-waller-factor=1.006)
        ```

        '''

        if isinstance(atom, CIFAtom):
            cell = atom.cell
            fractional_position = atom.fractional_position
            dwf = atom.dwf
            site_label = atom.site_label
            occupancy = atom.occupancy

        if isinstance(cell, Cell):
            self.cell = cell
        else:
            cell = np.array(cell, dtype=float)
            if cell.shape == (3, 3):
                self.cell = Cell(cell)
            elif cell.shape == (6,):
                dummy_cell = Cell(np.eye(3))
                self.cell = dummy_cell.new(cell)
            else:
                raise ValueError(
                    'Cell must be size (3,3) or (6,), got cell {cell!r} of shape {cell.shape}'.format(cell=cell))
        self.fractional_position = np.array(fractional_position, dtype=float)
        self.dwf = float(dwf)
        try:
            self.site_label = int(site_label)
        except ValueError:
            if isinstance(site_label, str):
                warn('String site labels may not be supported by MULTEM')
                self.site_label = site_label
            else:
                raise TypeError(
                    'Type of {site_label} (type {t}) is not supported for site_labels'.format(site_label=site_label,
                                                                                              t=type(site_label)))
        self.occupancy = float(occupancy)

        position = self.cell.cartesian_positions(self.fractional_position)
        super().__init__(symbol=symbol, charge=charge, position=position)

    def __str__(self):
        table = [[self.site_label, self.symbol, self.number, self.fractional_position[0], self.fractional_position[1],
                  self.fractional_position[2], self.dwf, self.occupancy, self.charge]]
        return '{rep}:\n{table}'.format(rep=repr(self), table=tabulate(table,
                                                                       headers=['Label', 'Symbol', 'Number', 'x/a',
                                                                                'y/b', 'z/c', 'B [Å]', 'Occupancy',
                                                                                'Charge']))

    @property
    def rms3d(self):
        return self.dwf / (8.0 * np.pi ** 2)

    @property
    def x(self):
        return self.position[0]

    @property
    def y(self):
        return self.position[1]

    @property
    def z(self):
        return self.position[2]

    def update_position(self):
        self.position = self.cell.cartesian_positions(self.fractional_position)


class Crystal(object):

    def __init__(self, atoms, cellpar, na=1, nb=1, nc=1):
        """
        Create a crystal object.

        A crystal consists of a unit cell and atoms with fractional coordinates. The atoms cartesian coordinates are calculated based on the unit cell (defined by the cell parameters `a`, `b`, `c`, `alpha`, `beta`, and `gamma`) on the fly.

        Parameters
        ----------
        atoms: list
            A list of CIFAtom objects
        cellpar: array-like
            A 6-length array
        na, nb, nc: int
            Bookkeeping numbers for keeping track of how many times the crystal has been "multiplied" - no real use except for when writing MULTEM files
        """
        self.atoms = atoms
        self.cellpar = np.array(cellpar, dtype=float)

        dummy_cell = Cell(np.eye(3))
        self.cell = dummy_cell.new(cell=self.cellpar)
        self.na = na
        self.nb = nb
        self.nc = nc

        site_labels = set([atom.site_label for atom in self.atoms])
        site_keys = site_labels
        self.site_dict = {}
        for label, key in zip(site_labels, site_keys):
            self.site_dict[label] = key

    def __repr__(self):
        if len(self.atoms) > 5:
            return '{self.__class__.__name__}(atom-list, {self.cellpar})'.format(self=self)
        else:
            return '{self.__class__.__name__}({self.atoms!r}, {self.cellpar})'.format(self=self)

    def __str__(self):
        cell_table = tabulate([self.cellpar],
                              headers=['a [Å]', 'b [Å]', 'c [Å]', 'alpha [deg]', 'beta [deg]', 'gamma [deg]'])
        orientation_table = tabulate([[name, x, y, z] for name, (x, y, z) in zip(['a', 'b', 'c'], self.cell.array)],
                                     headers=['', 'X [Å]', 'Y [Å]', 'Z [Å]'])
        atom_table = tabulate(
            [[atom.site_label, atom.symbol, atom.x, atom.y, atom.z, atom.dwf, atom.rms3d, atom.occupancy, atom.charge]
             for atom in self.atoms],
            headers=['Label', 'Symbol', 'x [Å]', 'y [Å]', 'z [Å]', 'B [Å]', 'RMS3D [Å]', 'Occ', 'Charge'])
        return 'Crystal with\n\nCell parameters:\n{cell_table}\n\nOrientation:\n{orientation_table}\n\nAtoms:\n{atom_table}'.format(
            cell_table=cell_table, orientation_table=orientation_table, atom_table=atom_table)

    def __mul__(self, other):
        other = np.array(other)
        new_cell = self.cell.new(
            [self.a * other[0], self.b * other[1], self.c * other[2], self.alpha, self.beta, self.gamma])
        new_atoms = []
        for atom in self.atoms:
            for i in range(other[0]):
                for j in range(other[1]):
                    for k in range(other[2]):
                        position = atom.position + np.array([i, j, k]) * atom.cell.lengths()
                        new_atom = CIFAtom(new_cell,
                                           fractional_position=new_cell.scaled_positions(position),
                                           symbol=atom.symbol,
                                           charge=atom.charge,
                                           dwf=atom.dwf,
                                           occupancy=atom.occupancy,
                                           site_label=atom.site_label)
                        new_atoms.append(new_atom)

        return Crystal(new_atoms, new_cell.cellpar(), na=self.na * other[0], nb=self.nb * other[1],
                       nc=self.nc * other[2])

    @property
    def a(self):
        return self.cellpar[0]

    @property
    def b(self):
        return self.cellpar[1]

    @property
    def c(self):
        return self.cellpar[2]

    @property
    def alpha(self):
        return self.cellpar[3]

    @property
    def beta(self):
        return self.cellpar[4]

    @property
    def gamma(self):
        return self.cellpar[5]

    @property
    def lx(self):
        return np.linalg.norm(self.cell.array, axis=0)[0]

    @property
    def ly(self):
        return np.linalg.norm(self.cell.array, axis=0)[1]

    @property
    def lz(self):
        return np.linalg.norm(self.cell.array, axis=0)[2]

    def min_inter_atomic_distance(self, axis=None, threshold=0):
        """
        Return the minimum interatomic distance above a threshold

        Parameters
        ----------
        axis: None or int, optional.
            The axis to consider (x=0, y=1, or z=2). If None (default), all axes are considered
        threshold: float
            Only consider interatomic distances above or at this threshold (default 0).

        Returns
        -------
        min_distance: float
            The minimum interatomic distance
        """
        min_distance = np.inf
        if axis is not None:
            if axis not in [0, 1, 2]:
                raise ValueError('Axis must be, 0, 1, or 2')
            for atom1 in self.atoms:
                for atom2 in self.atoms:
                    if atom1 is not atom2:
                        distance = abs(atom2.position[axis] - atom1.position[axis])
                        if threshold <= distance < min_distance: min_distance = distance
        else:
            for atom1 in self.atoms:
                for atom2 in self.atoms:
                    distance = abs(atom2.position - atom1.position)
                    if distance != 0 and distance < min_distance: min_distance = distance
        return min_distance

    def max_inter_atomic_distance(self, axis=None, threshold=None):
        """
        Return the maximum interatomic distance below a threshold

        Parameters
        ----------
        axis: None or int, optional.
            The axis to consider (x=0, y=1, or z=2). If None (default), all axes are considered
        threshold: float
            Only consider interatomic distances below or at this threshold (default is the largest cell dimension).

        Returns
        -------
        max_distance: float
            The maximum interatomic distance
        """
        if threshold is None:
            threshold = max([self.lx, self.ly, self.lz])
        max_distance = 0
        if axis is not None:
            if axis not in [0, 1, 2]:
                raise ValueError('Axis must be, 0, 1, or 2')
            for atom1 in self.atoms:
                for atom2 in self.atoms:
                    if atom1 is not atom2:
                        distance = abs(atom2.position[axis] - atom1.position[axis])
                        if distance != 0 and threshold >= distance > max_distance: max_distance = distance
        else:
            for atom1 in self.atoms:
                for atom2 in self.atoms:
                    distance = abs(atom2.position - atom1.position)
                    if distance != 0 and distance > max_distance: max_distance = distance
        return max_distance

    def add_atom(self, atom):
        """
        Append an atom to the crystal atom list

        Parameters
        ----------
        atom: CIFAtom
            The atom to be added

        """
        self.atoms.append(atom)

    def intifylabels(self):
        """
        Change the labels to integers.

        Changes the site labels of atoms to integers. The "key" is stored as well.

        Returns
        -------
        site_dict: dict
            A copy of the crystal site dictionary with the old site labels keys and the new labels as values
        """

        site_labels = set([atom.site_label for atom in self.atoms])
        site_keys = np.arange(0, len(site_labels), 1, dtype=int)
        site_dict = {}
        for label, key in zip(site_labels, site_keys):
            site_dict[label] = key

        for atom in self.atoms:
            atom.site_label = site_dict[atom.site_label]

        self.site_dict = site_dict
        return dict(self.site_dict)

    def write2mat(self, filename, dz, na=None, nb=None, nc=None):
        """
        Write the crystal to a .mat file for MULTEM

        Parameters
        ----------
        filename: str or Path
            The full path to save the .mat file
        dz: float
            The slice thickness (used for multislice simulations in MULTEM).
        na, nb, nc: int or None, optional.
            The number of replicated unit cells in a, b, and c directions (crystal directions). If None (default), the values stored in  the crystal object will be used.
        """

        filename = Path(filename)
        if any([isinstance(atom.site_label, str) for atom in self.atoms]):
            warn('String site labels detected for atoms in {self!r}, this might not be supported by MULTEM'.format(
                self=self))
        data = [[[atom.number], [atom.x], [atom.y], [atom.z], [atom.rms3d], [atom.occupancy], [atom.site_label],
                 [atom.charge]] if isinstance(atom.site_label, int) else [[atom.number], [atom.x], [atom.y], [atom.z],
                                                                          [atom.rms3d], [atom.occupancy],
                                                                          [[atom.site_label.encode('utf-8')]],
                                                                          [atom.charge]] for atom in self.atoms]

        matlabdict = {
            "spec_atoms": data,
            "spec_lx": self.lx,
            "spec_ly": self.ly,
            "spec_lz": self.lz,
            "spec_dz": dz,
            "a": self.a,
            "b": self.b,
            "c": self.c,
            "na": self.na,
            "nb": self.nb,
            "nc": self.nc,
            "site_keys": self.site_dict
        }
        savemat(filename.with_suffix('.mat'), matlabdict)
        print('Saved crystal to "{}"'.format(filename.with_suffix('.mat')))


class CIFfile(object):
    """
    A CIFfile converter object
    """

    def __init__(self, filename):
        """
        Create a CIFfile converter

        Parameters
        ----------
        filename: str or Path
            The path to the CIF file
        """
        self.filename = Path(filename)
        self.atoms = None
        self.crystal = None

    def read(self):
        """
        Read data from the CIF file
        """
        with ReadCif(str(self.filename)) as ciffile:
            for block in ciffile:
                # print(block)
                spacegroup = Spacegroup(int(block['_space_group_IT_number']))

                cellpar = np.array([
                    block['_cell_length_a'],
                    block['_cell_length_b'],
                    block['_cell_length_c'],
                    block['_cell_angle_alpha'],
                    block['_cell_angle_beta'],
                    block['_cell_angle_gamma'],
                ], dtype=float)

                try:
                    site_labels = block['_atom_site_label']
                except KeyError:
                    print('Could not get site labels from cif file')
                try:
                    occupancies = np.array(block['_atom_site_occupancy'], dtype=float)
                except KeyError:
                    print('Could not get occupancies from cif file')
                try:
                    fract_x = np.array(block['_atom_site_fract_x'], dtype=float)
                    fract_y = np.array(block['_atom_site_fract_y'], dtype=float)
                    fract_z = np.array(block['_atom_site_fract_z'], dtype=float)
                except KeyError:
                    warn('Could not get fractional coordinates from cif file, getting absolute coordinates instead.')
                    try:
                        x = np.array(block['_atom_site_cartn_x'], dtype=float)
                        y = np.array(block['_atom_site_cartn_y'], dtype=float)
                        z = np.array(block['_atom_site_cartn_z'], dtype=float)
                    except KeyError:
                        warn('Could not get absolute coordinates from cif file')
                        x = [np.nan] * len(block)
                        y = [np.nan] * len(block)
                        z = [np.nan] * len(block)
                    finally:
                        fract_x = x / cellpar[0]
                        fract_y = y / cellpar[1]
                        fract_z = z / cellpar[2]
                else:
                    x = fract_x * cellpar[0]
                    y = fract_y * cellpar[1]
                    z = fract_y * cellpar[2]
                finally:
                    x = x.T
                    y = y.T
                    z = z.T
                    fract_x = fract_x.T
                    fract_y = fract_y.T
                    fract_z = fract_z.T

                    positions = np.array([x, y, z])
                    fractional_positions = np.array([fract_x, fract_y, fract_z])
                try:
                    symbols = block['_atom_site_type_symbol']
                except KeyError:
                    print('Could not get atom site chemical symbols from cif file')
                try:
                    dwf = block['_atom_site_B_iso_or_equiv']
                except KeyError:
                    print('Could not get Debye-Waller factors from cif file')

                basis_atoms = []
                for label, occ, fx, fy, fz, symbol, B in zip(site_labels, occupancies, fractional_positions[0],
                                                             fractional_positions[1], fractional_positions[2], symbols,
                                                             dwf):
                    atom = CIFAtom(cellpar, symbol=symbol, occupancy=occ, fractional_position=(fx, fy, fz), dwf=B,
                                   site_label=label)
                    basis_atoms.append(atom)

                atoms = []
                for atom in basis_atoms:
                    equivalent_sites, kinds = spacegroup.equivalent_sites(atom.fractional_position, onduplicates='warn',
                                                                          occupancies=atom.occupancy)
                    for site in equivalent_sites:
                        position = site * cellpar[:3]
                        equivalent_atom = CIFAtom(cellpar, fractional_position=site, site_label=atom.site_label,
                                                  symbol=atom.symbol, dwf=atom.dwf, occupancy=atom.occupancy)
                        atoms.append(equivalent_atom)
                self.atoms = atoms
                self.crystal = Crystal(atoms, cellpar)


def convert_cif2mat(filename, dz, na=1, nb=1, nc=1, intify_labels=False):
    """
    Convert the crystal of a CIF file to .mat format

    Parameters
    ----------
    filename: str, Path
        Path to CIF file
    dz: float
        Slice thickness
    na, nb, nc: int, optional
        Number of crystal repetitions in crystal directions `a`, `b`, and `c`. Default is 1.
    intify_labels: bool, optional
        Whether to convert site labels into integers before writing the crystal to file or not. Default is False

    Returns
    -------
    ciffile: CIFfile
        The CIFfile with the atom and crystal data.
    """
    ciffile = CIFfile(filename)
    ciffile.read()
    if intify_labels:
        ciffile.crystal.intifylabels()
    ciffile.crystal = ciffile.crystal * [na, nb, nc]
    ciffile.crystal.write2mat(dz)
    return ciffile
