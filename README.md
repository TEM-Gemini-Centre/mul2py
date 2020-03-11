# mul2python

This is a python package for converting [MULTEM](https://github.com/Ivanlh20/MULTEM) results to [HyperSpy](https://github.com/hyperspy/hyperspy) signals and for subsequent plotting/visualisation. It requires specific fields to be saved to MATLAB structs from MULTEM, see more below.

## Disclaimer
This package is not meant for general use and is not tested. In particular, this package requires results from MULTEM in a specific format, and some knowledge of both MATLAB and python is needed to use it successfully. While the package is meant to be relatively simple to use, there are many pitfalls and details that must be considered, and care should be taken when using this package to generate results for scientific publication. This is especially true when it comes to overlaying atom positions on images and calibrating the signals from the simulation parameters. These things are only briefly tested and no guarantees are made. Possible problems in this regard involves trasposing (or not) the images relative to the atomic positions, reversial of scan directions, and image origins. Help is greatly appreciated, especially if unexpected results/errors are encountered.

## Installation
Install this package by downloading it and using pip. Navigate to where this file is downloaded and run pip in editable mode on the current directory:
```bash
$ cd <path to this directory>
$ pip install --editable .
```

## Workings
This package is divided in two sub-packages, `buildtools` and `exporttools`, for building and exporting results, respectively. The following lines will convert a MULTEM results file into a HyperSpy signal, and export a series of images that may be used to create videos:

### Pepare data
```Python
import matplotlib.pyplot as plt
import mul2py as m2p

#load signal
signal = m2p.build_ewrs('EWRS_results.mat')

#Add some metadata from original metadata
signal.metadata.add_dictionary({
    'SimulationParameters': {
        'E0': signal.original_metadata.Simulation_parameters.E0,
        'cond_lens_c_10': signal.original_metadata.Simulation_parameters.cond_lens_c_10,
        'cond_lens_c_30': signal.original_metadata.Simulation_parameters.cond_lens_c_30,
        'cond_lens_outer_aper_ang': signal.original_metadata.Simulation_parameters.cond_lens_outer_aper_ang,
        'nx': signal.original_metadata.Simulation_parameters.nx,
        'ny': signal.original_metadata.Simulation_parameters.ny,
        'spec_lx': signal.original_metadata.Simulation_parameters.spec_lx,
        'spec_ly': signal.original_metadata.Simulation_parameters.spec_ly,
        'spec_lz': signal.original_metadata.Simulation_parameters.spec_lz,
        'spec_dz': signal.original_metadata.Simulation_parameters.spec_dz,
        'thick': signal.original_metadata.Simulation_parameters.thick,
        'thick_type': signal.original_metadata.Simulation_parameters.thick_type,
        'spec_atoms': signal.original_metadata.Simulation_parameters.spec_atoms
    }
})

fig, ax = m2p.exporttools.make_image(signal, [0, 0, -3], False, 'annotate', 'mark_atoms',
                                     figure={'figsize': (4, 4),
                                             'dpi': 300})

#Make several images, convenient for making videos
counter = 0
thicknesses = [0, 9, 19, 29, 39] #Thicknesses to use
for s in signal:
    if signal.axes_manager.indices[-1] in thicknesses:
        fig, ax = m2p.exporttools.make_image(signal, None, True, 'annotate', 'mark_atoms',
                                             figure={'figsize': (4, 4),
                                                     'dpi': 300})
        
        frame_no = '{counter:0{pad:.0f}.0f}'.format(counter=counter, pad=3)
        fig.savefig('EWRS_{size:.0f}in_{dpi:.0f}dpi_{number}.png'.format(number=frame_no, size = fig.get_size_inches()[0], dpi = fig.dpi))
        
        plt.close(fig)
        
        counter += 1
    else:
        pass
```

### Commandline conversion
Alternatively, the results can be converted using the `convert_results.py` script:
```bash
$ source <path-to-suitabel-env>
$ python convert_results.py <path_to_data> <simulation_type>
```
`convert_results.py` will also set the `signal.metadata.SimulationParameters` similar to how it was done in the python cell above. 
## MULTEM
In MULTEM, the results should be saved in the following format
```MATLAB
input_multislice = JEM2100F_EWRS_setup("Al001_10x10x20.mat", 27.42); %Use premade function to load simulation parameters
original_input = input_multislice

centre_x = original_input.spec_lx/2;
centre_y = original_input.spec_ly/2;

step_x = original_input.spec_cryst_a/4;
step_y = original_input.spec_cryst_b/4;

width = 4; %Steps
height = 4; %Steps

xs = (centre_x - width / 2 * step_x : step_x : centre_x + width / 2 * step_x);
ys = (centre_y - height / 2 * step_y : step_y : centre_y + height / 2 * step_y);


results.input = original_input; %The input parameters used in the simulation
results.xs = xs; %The x-positions used for scanning - if applicable
results.ys = ys; %The y-positions used for scanning - if applicable
results.images = zeros(size(xs, 2), size(ys, 2), size(input_multislice.thick, 2), input_multislice.nx, input_multislice.ny); %Preallocate memory for image data
results.thicknesses = {}; %Cell for storing the thicknesses used at each x-y position.

%Loop over x-y- positions
counter = 1;
for i = 1:size(results.images, 1)
    for j = 1:size(results.images, 2)
        x = xs(i);
        y = ys(j);
        %shift beam
        input_multislice = original_input;
        input_multislice.iw_x = x;
        input_multislice.iw_y = y;
        
        file_name = sprintf("%s_%i_%i", result_name, i, j); %Name for storing individual output file if an exception is thrown.

        %Run simulation
        clear il_MULTEM;
        fprintf('Simulating exit wave (real-space) stack %i of %i: (x,y) = (%f,%f)\r', counter, size(xs, 2)*size(ys, 2), input_multislice.iw_x, input_multislice.iw_y);
        
        tic;
        output_multislice = il_MULTEM(system_conf, input_multislice);
        toc;
        results.thicknesses{i, j} = output_multislice.thick;

        try
            for t = 1:size(output_multislice.data, 2)
                results.images(i, j, t, :, :) = output_multislice.data(t).m2psi_tot;
            end
        catch ME
            fprintf('Exception for i=%i, j=%i, and t=%i. Data size: (%s)', i, j, t, strip(sprintf('%i,', size(output_multislice.data)), 'right', ','));
            save(sprintf("%s/%s_output.mat", results_path, file_name), 'output_multislice', '-v7.3');
            save(sprintf("%s/%s_results.mat", results_path, result_name), 'results', '-v7.3');
            rethrow(ME)
        end        
        counter=counter+1;
    end
end
results.dx = output_multislice.dx;
results.dy = output_multislice.dy;

save(sprintf("%s/%s_results.mat", results_path, result_name), 'results', '-v7.3');
```
