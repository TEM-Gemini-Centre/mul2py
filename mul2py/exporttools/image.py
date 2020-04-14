import matplotlib.pyplot as plt
import numpy as np

def make_image(signal, inav=None, transpose=False, mark_atoms=True, annotations=True, **kwargs):
    """Make an image of simulation results.

    Shows a frame of the signal and adds annotations and markings to it. The frame can be selected from a provided 'inav' index, or by picking the current axes_manager indices.

    Parameters
    ----------
    signal : hyperspy.api.signals.Signal2D
        The signal to show an image from.
    inav : None or array-like, optional.
        The frame from signal to show. If None, then the current axes_manager.indices is picked.
    transpose : bool, optional.
        Whether to transpose the data before plotting. Default is False.
    mark_atoms : bool, optional.
        Whether to mark atoms above the current image slice in the model or not. The appearance of the markings is controlled by kwargs['markers']. Default is True.
    annotations : bool, optional.
        Whether to annotate the image according to the kwargs['annotate'] parameters. Default is True

    Other Parameters
    ----------------
    **kwargs : optional keyword arguments passed on to sub-calls, optional.
        kwargs['figure'] is passed to plt.figure(). Default is {'dpi': 72, 'figsize': (2,2), 'frameon': False
        kwargs['ax'] is passed to fig.add_axes(). Default is {'position': [0,0,1,1], 'xticks': [], 'yticks': [], 'frameon': False}
        kwargs['imshow'] is passed to plt.imshow(). Default is {'extent': (min(y), max(y), min(x), max(x)), 'vmin': min(data), 'vmax':max(data)}
        kwargs['annotate'] is passed to plt.annotate(). Default is {'color': 'w', 'ha':'left', 'va':'top', 'xy':(0.02, 0.98), 'xycoords': 'axes fraction', 'format_spec':'.2f', 'text':None}. 'text':None will annotate the image with corresponding navigation space coordinates.
        kwargs['markers'] is passed to plt.scatter(). Default is {'edgecolors': Z-value of atoms, 's':10, 'lw':0.4, 'marker':'o', 'alpha':0.5, 'facecolors': 'none'}.

    Returns
    -------
    fig, ax : the figure and axis of the image.
    """
    # Set up figure
    fig_kwargs = {
        'dpi': 72,
        'figsize': (2, 2),
        'frameon': False
    }
    try:
        fig_kwargs.update(kwargs.pop('figure'))
    except KeyError:
        pass
    finally:
        fig = plt.figure(**fig_kwargs)

    # Set up axis
    ax_kwargs = {
        'position': [0, 0, 1, 1],
        'xticks': [],
        'yticks': [],
        'frameon': False
    }
    try:
        ax_kwargs.update(kwargs.pop('ax'))
    except KeyError:
        pass
    finally:
        ax_pos = ax_kwargs.pop('position')
        ax = fig.add_axes(ax_pos, **ax_kwargs)

    # Set up navigation indices
    if inav is not None:
        inav = [ax.high_index + i if i < 0 else i for ax, i in
                zip(signal.axes_manager.navigation_axes, inav)]  # Convert negative indices to high indices.
        signal.axes_manager.trait_set(indices=inav)
    inav = signal.axes_manager.trait_get('indices')['indices']

    # Make sure the signal is an image and set up image extent
    sig_axes = [ax.name for ax in signal.axes_manager.signal_axes]
    assert len(sig_axes) == 2, 'There are {} signal axes, expected 2!'.format(len(sig_axes))
    xs = signal.axes_manager[sig_axes[0]].axis
    ys = signal.axes_manager[sig_axes[1]].axis
    extent = (np.min(xs), np.max(xs), np.min(ys), np.max(ys))

    # Set up image plot
    imshow_kwargs = {
        'extent': extent,
        'vmin': np.min(signal.inav[inav].data),
        'vmax': np.max(signal.inav[inav].data),
        'origin': 'lower',
    }
    try:
        imshow_kwargs.update(kwargs.pop('imshow'))
    except KeyError:
        pass
    finally:
        if transpose:
            ax.imshow(signal.inav[inav].data.T, **imshow_kwargs)
        else:
            ax.imshow(signal.inav[inav].data, **imshow_kwargs)

    # Set up frame annotation
    if annotations:
        annotate_kwargs = {
            'color': 'w',
            'ha': 'left',
            'va': 'top',
            'xy': (0.02, 0.98),
            'format_spec': '.2f',
            'xycoords': 'axes fraction',
            'text': None
        }
        try:
            annotate_kwargs.update(kwargs.pop('annotate'))
        except KeyError:
            pass
        finally:
            if annotate_kwargs['text'] is None:
                format_spec = annotate_kwargs.pop('format_spec')
                annotate_kwargs.update({'text': ', '.join(['{ax.name}={val:{format_spec}} {ax.units}'.format(ax=ax,
                                                                                                             val=ax.scale * ax.index + ax.offset,
                                                                                                             format_spec=format_spec)
                                                           for ax in signal.axes_manager.navigation_axes])})

            text = annotate_kwargs.pop('text')
            xy = annotate_kwargs.pop('xy')
            ax.annotate(text, xy, **annotate_kwargs)

    # Set up atom markings
    if mark_atoms:
        atoms = np.array(signal.metadata.SimulationParameters.spec_atoms).T
        atoms = atoms[
            atoms[:, 3] < signal.axes_manager['z'].offset + signal.axes_manager['z'].index * signal.axes_manager[
                'z'].scale]  # Only mark atoms "above" the current slice.

        atoms = set([(atom[0], atom[1], atom[2]) for atom in np.array(
            signal.metadata.SimulationParameters.spec_atoms.T)])  # Get only the "unique" atoms in the x-y plane
        atoms = np.array([np.array(atom) for atom in atoms])  # Change into numpyarray

        marker_kwargs = {
            's': 10,
            'lw': 0.4,
            'marker': 'o',
            'alpha': 0.5,
            'edgecolors': None,
            'cmap': 'RdBu',
            'facecolors': 'none'
        }
        try:
            marker_kwargs.update(kwargs.pop('markers'))
        except KeyError:
            marker_kwargs
        finally:
            if marker_kwargs['edgecolors'] is None:
                cmap = plt.cm.get_cmap(marker_kwargs.pop('cmap'))
                marker_kwargs['edgecolors'] = cmap(atoms[:, 0])
            marker_kwargs['edgecolors']
            if marker_kwargs['facecolors'] == 'edgecolors':
                marker_kwargs['facecolors'] = marker_kwargs['edgecolors']
            ax.scatter(atoms[:, 1], atoms[:, 2], **marker_kwargs)

    return fig, ax

