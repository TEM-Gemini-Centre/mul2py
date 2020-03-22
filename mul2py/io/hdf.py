import h5py as hdf
import numpy as np
import tabulate as tab


class Error(Exception):
    pass


class HDFContentError(Error):
    pass


class HDFReaderError(Error):
    pass


class HDFContent(object):
    """Object for storing hdf content"""

    def __init__(self, name, content, hdf_file):
        self._file = hdf_file
        self._name = name
        self._content = content
        if isinstance(self._content, hdf.Reference):
            self._is_reference = True
        else:
            self._is_reference = False
            self._depth = len(self._content)

        if isinstance(self._content, hdf.Group):  # If the content is a group, unpack it
            for subcontent in self._content:
                self[subcontent] = HDFContent(subcontent, self._content[subcontent], self._file)

        elif isinstance(self._content,
                        hdf.Dataset):  # if the content is a dataset, do nothing unless the dataset contains general
            # objects (i.e. references or groups)
            self._shape = self._content.shape
            self._content = np.array(self._content)
            if self._content.dtype == 'O':  # if the dataset is "objects", unpack further
                for i, data in enumerate(self._content.flat):
                    attribute_name = '{name}{i}'.format(name=self._name, i=i)
                    self[attribute_name] = HDFContent('resolved_reference', data, self._file)
        elif isinstance(self._content, hdf.Reference):
            self[self._name] = HDFContent(self._name, self._file[self._content], self._file)

    @property
    def content(self):
        return self._content

    @content.setter
    def content(self, content):
        self._content = content

    @property
    def name(self):
        return self._name

    @property
    def file(self):
        return self._file

    @property
    def is_reference(self):
        return self._is_reference

    @property
    def shape(self):
        return self._shape

    @property
    def depth(self):
        return self._depth

    def __repr__(self):
        return '{self.__class__.__name__}({self._name!r}, {self._content!r})'.format(self=self)

    def __str__(self):
        table = tab.tabulate([[attribute.name, str(attribute.content)] for attribute in self],
                             headers=['name', 'content'])
        return 'HDF content {self.name}:\n{table}'.format(self=self, table=table)

    def __getitem__(self, item):
        if isinstance(item, self.__class__):
            item = item.name
        return getattr(self, item)

    def __setitem__(self, key, value):
        if isinstance(key, self.__class__):
            key = self.name
        if key[0] == '_':
            raise KeyError('Attempting to set protected value "{}"'.format(key))
        else:
            setattr(self, key, value)

    def __iter__(self):
        for key in self.__dict__.keys():
            if any([key[0] == invalid_prefix for invalid_prefix in ('_', '#')]):
                pass
            else:
                yield self[key]

    def __call__(self):
        try:
            return self.content
        except HDFContentError:
            return self

    def __len__(self):
        return len(list(iter(self)))

    def remove_excess_dimensions(self):
        """Removes excess dimensions of content"""
        if len(self) > 0:
            for attribute in self:
                # print('Removing excess dimensions of {self.name}.{attribute.name}'.format(self=self, attribute=attribute))
                if isinstance(attribute, self.__class__):
                    # print('Attribute is type {self.__class__.__name__}'.format(self=self))
                    try:
                        attribute.remove_excess_dimensions()
                    except HDFContentError as e:
                        # print('Excess dimensions of {att} could not be removed'.format(att=attribute))
                        raise e
        else:
            if self.shape == (1, 1):  # 0D array
                # print('0D array')
                self.content = self.content.ravel()[0]  # 2D -> 0D
            elif len(self.shape) == 2 and self.shape[0] > 1 and self.shape[1] == 1:  # 1D array
                # print('1D array')
                self.content = self.content.ravel()  # 2D -> 1D
            elif len(self.shape) == 2 and self.shape[0] > 1 and self.shape[1] > 1:  # 2D array
                # print('2D array')
                pass  # Leave 2D arrays as is

    def content2dict(self):
        """Return contents as a dictionary"""
        if len(self) > 0:
            dictionary = {}
            for attribute in self:
                dictionary[attribute.name] = attribute.content2dict()
            return dictionary
        else:
            return self.content


class HDFReader(object):
    """Object for reading hdf files"""

    def __init__(self, hdf_file):
        if isinstance(hdf_file, hdf.File):
            self._file = hdf_file
        else:
            self._file = hdf.File(hdf_file, 'r')  # Try to treat the argument as a string/path/etc

        self.read()

    def __repr__(self):
        return '{self.__class__.__name__}({self.file})'.format(self=self)

    def __str__(self):
        table = tab.tabulate([[attribute.name, str(attribute)] for attribute in self], headers=['name', 'content'])
        return 'HDF file with content:\n{table}'.format(table=table)

    def __getitem__(self, item):
        if isinstance(item, HDFContent):
            item = item.name
        return getattr(self, item)

    def __setitem__(self, key, value):
        if isinstance(key, HDFContent):
            key = key.name
        if key[0] == '_':
            raise KeyError('Attempting to set protected value "{}"'.format(key))
        else:
            setattr(self, key, value)

    def __iter__(self):
        for key in self.__dict__.keys():
            if any([key[0] == invalid_prefix for invalid_prefix in ('_', '#')]):
                pass
            else:
                yield self[key]

    @property
    def file(self):
        return self._file

    def read(self):
        for field in self.file:
            self[field] = HDFContent(field, self.file[field], self.file)

    def close(self):
        self.file.close()

    def remove_excess_dimensions(self, *args):
        """Removes excess dimensions of the data

        Other Parameters
        ----------------
        *args: str
            The contents to remove the excess dimensions of

        """
        if len(args) > 0:  # remove excess dimensions of attributes named "arg" in args
            for arg in args:
                if isinstance(arg, HDFContent):
                    # print('Removing excess dimensions of content {item.name}'.format(item=self[arg]))
                    self[arg].remove_excess_dimensions()
        else:  # remove excess dimensions of all attributes
            for attribute in self:
                self.remove_excess_dimensions(attribute)

    def content2dict(self):
        dictionary = {}
        for attribute in self:
            dictionary[attribute.name] = attribute.content2dict()
        return dictionary
