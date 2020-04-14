import pathlib
import numpy as np
import matplotlib.pyplot as plt
from .image import make_image as _make_image


def make_movie(filepath, signal, inavs, annotate=True, show_progress=False, mark_atoms=True, **kwargs):
    """Export images from signal for movie creation

    Parameters
    ---------
    filepath : str
        the output file path
    signal : HyperSpy signal
        the signal
    inavs : list of N lists,
        The navigation axis indexes to export data from. Should be a list of N lists, where N is the signal navigation dimension.
    annotate : bool, optional
        Whether to annotate the images with navigation axis details. Default is True.
    show_progress : bool, optional
        Whether to show the progress
    mark_atoms : bool, optional
        Whether to mark the atoms above each slice. Default is True.
    **kwargs: optional keyword arguments passed to MULTEM2python.exporttools.image.make_image()
    """

    filepath = pathlib.Path(filepath)
    number_of_frames = np.array([len(inav) for inav in inavs]).cumprod()[-1]  # The number of images to make

    counter = 0
    for inav in zip(inavs):
        if show_progress:
            print('Making image {} of {}'.format(counter + 1, number_of_frames))
        inav = inav[0]  # the zip is a tuple of one list - extract the list
        fig, ax = _make_image(signal, inav=inav, transpose=False, mark_atoms=mark_atoms, annotations=annotate, **kwargs)
        size = fig.get_size_inches()[0]
        dpi = fig.get_dpi()
        number = '{counter:0{pad:.0f}.0f}'.format(counter=1, pad=int((np.floor(np.log10(number_of_frames)) + 2)))
        fig.savefig(filepath.with_name(
            '{stem}_{size:.0f}in_{dpi:.0f}dpi_{number}'.format(stem=filepath.stem, number=number, size=size,
                                                               dpi=dpi)), dpi=dpi, figsize=(size, size))
        plt.close(fig)
        counter += 1
