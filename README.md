# mul2py

This is a python package for converting [MULTEM](https://github.com/Ivanlh20/MULTEM) results to [HyperSpy](https://github.com/hyperspy/hyperspy) signals and for subsequent plotting/visualisation. It requires specific fields to be saved to MATLAB structs from MULTEM, see more [here](#5-the-easy-mat-format).

## Content
  - [1 Disclaimer](#1-disclaimer)
  - [2 Installation](#2-installation)
  - [3 Workings](#3-workings)
    - [3.1 Creating signals](#31-creating-signals)
      - [3.1.1 Easy way (.ecmat)](#311-easy-way)
      - [3.1.2 Robust way (.mat)](#312-robust-way)
    - [3.2 Making annotated images and movies](#32-making-images-and-movies)
  - [4 Commandline conversion](#4-commandline-results-conversion)
  - [5 The easy-mat format](#5-the-easy-mat-format)
    - [5.1 Easy-mat example](#51-matlab-multem-example)
  - [6 Examples](#6-examples)
    - [6.1 Example model script](#61-making-a-model)
    - [6.2 Example STEM script](#62-stem-simulation-script)
    - [6.3 Example SLURM job](#63-slurm-job-script)

## 1 Disclaimer
This package is not meant for general use and is not tested. In particular, this package requires results from MULTEM in a specific format, and some knowledge of both MATLAB and python is needed to use it successfully. While the package is meant to be relatively simple to use, there are many pitfalls and details that must be considered, and care should be taken when using this package to generate results for scientific publication. This is especially true when it comes to overlaying atom positions on images and calibrating the signals from the simulation parameters. These things are only briefly tested and no guarantees are made. Possible problems in this regard involves trasposing (or not) the images relative to the atomic positions, reversial of scan directions, and image origins. Help is greatly appreciated, especially if unexpected results/errors are encountered.

[Up](#content)

## 2 Installation
Install this package by downloading it and using pip. Navigate to where this file is downloaded and run pip in editable mode on the current directory:
```bash
$ cd <path to this directory>
$ pip install --editable .
```
Note that the trailing `.` is very important.

[Up](#content)

## 3 Workings
This package is divided in several subpackages. The [user guide](https://github.com/TEM-Gemini-Centre/mul2py/tree/master/UserGuide.md) is a good place to start to learn about these packages, but it is far from an extensive documentation. In short, mul2py provides the subpackages `buildtools`. `exporttools`, and `io` that can be used for building, exporting, and reading/writing, respectively.

[Up](#content)

### 3.1 Creating signals
To build a hyperspy signal from MULTEM results, there are two main options, using `mul2py.buildtools.builders.make_signal()` or by using `mul2py.buildtools.builders.load_output()` together with `mul2py.buildtools.builders.load_input()` to generate useful metadata.

[Up](#content)

#### 3.1.1 Easy way
To use `mul2py.buildtools.builders.make_signal()`, the results file from MULTEM must be given in a ".ecmat" (pronounced easy-mat) file. This is simply a normal MATLAB v7.3 mat file with data stored in very specific fields (see more later). This allows calibrated signals with meaningful metadata to be generated by:
```Python
import mul2py as m2p
signal = m2p.buildtools.make_signal('my-results.ecmat')
``` 

[Up](#content)

#### 3.1.2 Robust way
If the data from MULTEM is instead given in 'my-output.mat' and 'my-input.mat'. that contains the fields/groups 'output_multislice' and 'input_multislice' at their roots, the signals and metadata can be generated by:
```Python
import mul2py as m2p
metadata = m2p.buildtools.load_input('my-input.mat')
signal, output = m2p.buildtools.load_output('my-output.mat')
```
These lines will produce an uncalibrated an raw image stack in the `signal` variable, while `metadata` will be a `dict` that contains all the input parameters to the simulation. `output` will be a `HDFContent` object that contains the fields of the raw file and can be used for debugging or further metadata or data tuning.

[Up](#content)

### 3.2 Making images and movies
Images with annotated atom positions can be made by
```Python
import mul2py as m2p
signal = m2p.buildtools.make_signal('my-results.ecmat')
fig, ax = m2p.exporttools.make_image(signal)
```
The appearance of the atomic positions can be tuned by providing a keyword dictionary that is passed to `matplotlib.pyplot.scatter()`:
```Python
import mul2py as m2p
signal = m2p.buildtools.make_signal('my-results.ecmat')
fig, ax = m2p.exporttools.make_image(signal, markers={'cmap': 'magma'})
```

In some cases, it is also useful to create several images through the image stacks to visualise how the iamges evolve through the thickness and for different beam positions. This can be done by
```Python
import mul2py as m2p
import numpy as np
signal = m2p.buildtools.make_signal('my-results.ecmat')
inavs = np.arange(0, len(signal)-1, 2) #Take every second slice. NB! Only for 3D simulations. For 4D or 5D simulations, you must provide a 2D or 3D list of navigation indices.
m2p.exporttools.make_movie('path-to-my-image-series', signal, inavs=inavs)
```

[Up](#content)

## 4 Commandline results conversion
Results can also be converted using the `convert_results.py` script from the commandline:
```bash
$ source <path-to-suitable-env>
$ python convert_results.py <path_to_data> <simulation_type>
```

[Up](#content)

## 5 The easy-mat format
In MULTEM, the results should be saved in the following format for `mul2py.buildtools.make_signal()` to work properly:
  - results
    - input: The input `struct` passed to the simulation
    - system: The system configuration `struct` passed to the simulation
    - xs: The x-positions of the beam
    - ys: The y-positions of the beam
    - images: A nD stack of images (the actual results)
    - thicknesses: A cell of thicknesses at each x-y-position
    - thick: The thicknesses from the last output (i.e. the output from `il_MULTEM()`)
    - dx: Image resolution in x-direction
    - dy: Image resolution in y-direction
    - elapsed_time: The elapsed time of the complete MATLAB work (from start of the script to the end right before storing data)

In MATLAB, this is done by creating a struct called `results` with the various fields:
```MATLAB
input_multislice; %The input parameter to il_MULTEM
system_conf; %The system configuration passed to il_MULTEM
output_multislice = il_MULTEM(system_conf, input_multislice);

results.input=input_multislice;
results.system = system_conf;
results.xs = xs;
results.ys = ys;
results.images = image_stack;
results.thicknesses = {};
results.thick = output_multislice.thick;
results.dx = output_multislice.dx;
results.dy = output_multislice.dy;
results.elapsed_time = seconds(stop-start);
```

The units are always given as Å, except for `dx` and `dy`. For real-space simualtion images such as EWRS, STEM, and HRTEM simulations, these values are given as Å. For diffraction-space simulations such as CBED, they have units 1/Å.

[Up](#content)

### 5.1 MATLAB MULTEM Example
```MATLAB
%%
clear all
clc
%% System configuration
gpu = true; %Use GPU?
if gpu
    system_conf.precision = 1;                                              % eP_Float = 1, eP_double = 2
    system_conf.device = 2;                                                 % eD_CPU = 1, eD_GPU = 2
    system_conf.cpu_nthread = 5; 			                                % Does the number of CPU threads matter when running on GPU? EXPERIMENT!!
    system_conf.gpu_device = 1;				                                % MULTEM can only use one GPU device at the time? Only ask for a single GPU from IDUN, and use this.
else
    system_conf.precision = 1;                                              % eP_Float = 1, eP_double = 2
    system_conf.device = 1;                                                 % eD_CPU = 1, eD_GPU = 2
    system_conf.cpu_nthread = 4; 
    system_conf.gpu_device = 0;
end

%% Timestamp
start_time = datetime('now','TimeZone','local');                            % Start tracking time
fprintf("Starting simulation script at %s\n", start_time);                  % Print the start time

%% output_details
simulation_name = "EWRS_test";                                              % Make a name for the simulation
output_path = ".";                                                          % Define an output path
mkdir(output_path);                                                         % Make the output directory

%% Make simulation parameters
input_multislice = EWRS_setup("test_model_L_10x10x20.mat", 15, "phonons", 1, "nx", 8, "ny", 16, "instrument", "2100F", "multem_path", "C:\Program Files\MULTEM\MULTEM_binary");
original_input = input_multislice;                                          % Backup the original input parameters
results.system = system_conf;                                               % Store the system configuration in the `results` structure

%% Set up scan pattern
centre_x = original_input.spec_lx/2;                                        % Scan centre
centre_y = original_input.spec_ly/2;                                        % Scan centre

scanning_width = original_input.spec_cryst_a;                               % Scan width in Å
scanning_height = original_input.spec_cryst_b;                              % Scan height in Å

scanning_ns_x = 2;                                                          % Number of scan points
scanning_ns_y = 3;                                                          % Number of scan points

x0 = centre_x - scanning_width/2;                                           % Lower left corner of scan (model viewed from bottom?)
y0 = centre_y - scanning_height/2;                                          % Lower left corner of scan (model viewed from bottom?)

xe = x0 + scanning_width;                                                   % Upper right corner of scan (model viewed from bottom?)
ye = x0 + scanning_height;                                                  % Upper right corner of scan (model viewed from bottom?)


xs = linspace(x0, xe, scanning_ns_x);                                       % Beam scan positions
ys = linspace(y0, ye, scanning_ns_y);                                       % Beam scan positions

%% Loop through x and y positions
results.input = original_input;                                             % Store input parameters
results.xs = xs;                                                            % Store beam positions
results.ys = ys;                                                            % Store beam positions
results.images = zeros(input_multislice.nx, input_multislice.ny, size(xs, 2), size(ys, 2), size(input_multislice.thick, 2)); %Predefine image stack
results.thicknesses = {};                                                   % Predefine thickness cell

%% Loop through beam positions and do a simulation at each position
counter = 1;
for i = 1:size(results.xs, 2)
    x = xs(i);
    for j = 1:size(results.ys, 2)
        y = ys(j);
        %shift beam
        input_multislice = original_input;
        input_multislice.iw_x = x;
        input_multislice.iw_y = y;
        
        %Run simulation
        fprintf("Simulating EWRS stack %i of %i: (x,y) = (%f,%f)\r", counter, length(xs) * length(ys), input_multislice.iw_x, input_multislice.iw_y);
        clear il_MULTEM;
        tic;
        output_multislice = il_MULTEM(system_conf, input_multislice); 
        toc;
        results.thicknesses{i, j} = output_multislice.thick;
        try
            for t = 1:length(output_multislice.data)
                results.images(:, :, i, j, t) = transpose(output_multislice.data(t).m2psi_tot); %Note that the images must be transposed to align properly in python (using different column-order)
            end
        catch ME
            fprintf("Exception for i=%i, j=%i, and t=%i. Data size: (%s)", i, j, t, strip(sprintf("%i,", size(output_multislice.data)), "right", ","));
            save(sprintf("%s/%s_%i_%i_%i_output.mat", output_path, simulation_name, i, j, t), "output_multislice", "-v7.3");
            save(sprintf("%s/%s_%i_%i_%i_results.ecmat", output_path, simulation_name, i, j, t), "results", "-v7.3");
            rethrow(ME)
        end        
        counter=counter+1;
    end
end
%% Set some additional details
results.thick = output_multislice.thick;
results.dx = output_multislice.dx;
results.dy = output_multislice.dy;

end_time = datetime('now','TimeZone','local');
fprintf("Simulation finished at %s\n", end_time);
results.elapsed_time = seconds(end_time - start_time);

%% Save the data
save(sprintf("%s/%s_results.ecmat", output_path, simulation_name), "results", "-v7.3");
```

[Up](#content)

## 6 Examples
The general workflow of MULTEM is
  1. Make a model in your favourite tool and convert it to a .mat file
        - Use e.g. `mul2py.io.models.save_multem_model()` to convert it.
  2. Write a MATLAB script that specifies simulation parameters and loads this model, called e.g. "my_simulation.m"
        - This can be done by using the setup-functions provided at [https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/matlab](https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/matlab)
  3. Run a SLURM job on the cluster
        - Run the simulation in MATLAB: `matlab -nodisplay -nodesktop -nosplash -r "my_simulation"`
        - Navigate to where the results (named e.g. "my_results.ecmat") were stored and convert your results using: `python convert_ecmat.py my_results.ecmat`
        - Download the resulting .hspy file and enjoy!

Several complete examples are also provided in [https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/examples](https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/examples).

### 6.1 Making a model
To make a model based on a crystal unit cell defined in "Al.cif", simply load the cif file with [ASE](https://wiki.fysik.dtu.dk/ase/) and replicate it as you want, then save it as a .mat file for MULTEM using `mul2py.io.models.save_multem_model()`
```Python
from ase.io import read
from mul2py.io.models import save_multem_model

na, nb, nc = 10, 10, 20
Al = read('Al.cif')
slab = Al*[na, nb, nc]

#Define debye-waller factors for the elements in the model in a dictionary. This is required by MULTEM
dwfs = {
    13: 0.7993,
}

#Define the slice thickness you want to slice the model with
a, b, c, alpha, beta, gamma = Al.get_cell_lengths_and_angles() #Get lattice parameters
dz = c/2 #Slice the model in a meaningful way (for fcc, a slice thickness of a/2 will produce a slicing where only one atomic layer is inside the slice)

save_multem_model('Al_{nx:0f}x{ny:0f}x{nz:0f}.mat'.format(nx=na, ny=nb, nz=nc), slab, B=dwfs) #Save the model.
```

### 6.2 STEM simulation script
As an example, a STEM simulation script "STEMcpu.m" may look like this:
```MATLAB
%%
clear all
clc

%% Set up paths
multem_path =  '/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM' %Path to MULTEM installation
addpath(char(sprintf("%s/crystalline_materials", multem_path))); % Add specific subfolders to the path
addpath(char(sprintf("%s/matlab_functions", multem_path))); % Add specific subfolders to the path
addpath(char(sprintf("%s/mex_bin", multem_path))); % Add specific subfolders to the path

%% System configuration
system_conf.precision = 1;                           % eP_Float = 1, eP_double = 2
system_conf.device = 1;                              % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 16; 
system_conf.gpu_device = 0;

%% Timestamp
start_time = datetime('now','TimeZone','local');
fprintf("Starting simulation script at %s\n", start_time);

%% Define detector collection semi-angles
clear collection_angles
collection_angles(1).inner_ang = 0; %semi-angle (mrad)
collection_angles(1).outer_ang = 40; %semi-angle (mrad)
collection_angles(2).inner_ang = 48; %semi-angle (mrad)
collection_angles(2).outer_ang = 200; %semi-angle (mrad)

%% Define beam convergence semi-angle
convergence_angle = 27; %semi-angle (mrad)

input_multislice = STEM_setup("test_model_L_10x10x20.mat", convergence_angle, collection_angles, "phonons", 1, "nx", 8, "ny", 16, "instrument", "ARM200F", "multem_path", multem_path);

%% Store input parameters in results struct
results.input = input_multislice;
results.system = system_conf;

%% Run simulation
clear il_MULTEM;
tic;
output_multislice = il_MULTEM(system_conf, input_multislice); 
toc;

%% Construct results.images
n_t = length(output_multislice.data); %The number of thicknesses
detectors = length(output_multislice.data(1).image_tot); % The number of detectors

results.images = zeros(input_multislice.scanning_ns, input_multislice.scanning_ns, n_t, detectors);

for t = 1:n_t
    for d = 1:detectors
        results.images(:, :, t, d) = transpose(output_multislice.data(t).image_tot(d).image);
    end
end

%% Set some additional details
results.thick = output_multislice.thick;
results.dx = output_multislice.dx;
results.dy = output_multislice.dy;

end_time = datetime('now','TimeZone','local');
fprintf("Simulation finished at %s\n", end_time);
results.elapsed_time = seconds(end_time - start_time);

%% Save the data
save("/lustre1/work/emilc/MULTEM/Test/STEMcpu_results.ecmat", "results", "-v7.3"); %Choose a sensible directory to save your results.
``` 
 
Example MATLAB scripts can be found at [https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/matlab/example_multem_scripts](https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/matlab/example_multem_scripts). These scripts are meant for testing mainly, but they show the general idea quite nicely as well.

### 6.3 SLURM job script
When running on the cluster, you must choose the correct queue that matches your `system_conf` variable in the matlab script. As an example, the following slurm code will run the above simulation and then convert it to a hyperspy signal.

```bash
#!/bin/bash

#SBATCH --partition=CPUQ
#SBATCH --time=01-20:0:00
#SBATCH --job-name="STEMcpu"
#SBATCH --output=STEMcpu-%A.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64000
#SBATCH --account=share-nv-fys-tem

echo "we are running from this directory: $SLURM_SUBMIT_DIR"
echo "The name of the job is: $SLURM_JOB_NAME"

echo "The job ID is $SLURM_JOB_ID"
echo "The job was run on these nodes: $SLURM_JOB_NODELIST"
echo "Number of nodes: $SLURM_JOB_NUM_NODES"
echo "We are using $SLURM_CPUS_ON_NODE cores"

echo "We are using $SLURM_CPUS_ON_NODE cores per node"
echo "Total of $SLURM_NTASKS cores"

module load foss/2016a #Necessary for CUDA
module load CUDA/8.0.61 #Necessary for MULTEM, even if no GPU is being used I think
module load MATLAB/2017a #Necessary for MULTEM

echo "Running STEM simulation"
matlab -nodisplay -nodesktop -nosplash -r "STEMcpu" #Run the script called `STEMcpu.m`

echo "Converting results to HyperSpy format using mul2py"
module load GCCcore/.8.2.0 Python/3.7.2 #Necessary for mul2py
source /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py-env/bin/activate #Activate environment with mul2py in project folder
python /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py/mul2py/examples/convert_ecmat.py /lustre1/work/emilc/MULTEM/Test/STEMcpu_results.ecmat #Run convert_ecmat.py script. The last path here must match the path you chose to store your results at in the matlab script.

scontrol show job ${SLURM_JOB_ID} -d #Print job details
```

