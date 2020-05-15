# mul2py User guide
This is a user guide to introduce mul2py and [MULTEM](https://github.com/Ivanlh20/MULTEM) simulations to new users. It is not intended as an extensive or complete guide, documentation, or reference. However, it outlines the more important aspects of the package and how it can be used together with MULTEM
## Content
  - [1 Introduction](#1-introduction)
  - [2 Installation](#2-installation)
    - [2.1 Clone mul2py](#21-step-1-clone-mul2py)
      - [2.1.1 Update mul2py](#211-to-update-mul2py)
    - [2.2 Virtual environments](#22-step-2-activatecreate-environment)
    - [2.3 Install mul2py](#23-step-3-add-mul2py-to-environment)
  - [3 Convert MULTEM results](#3-convert-multem-results)
    - [3.1 Working concept](#31-working-concept)
  - [4 Plotting images and making movies](#4-plotting-images-and-exporting-data-for-movies)
  - [5 Making models](#5-making-models)
  - [6 Running MULTEM simulations](#6-running-multem-simulations)
    - [6.1 Running simulations from shell](#61-running-simulations-from-shell)
    - [6.2 Calling matlab scripts from shell](#62-calling-matlab-script-from-shell)
        - [6.2.1 Running run-functions through a script](#621-running-run_simulation-type_simulation-through-a-script)
        - [6.2.2 Running simulations using set-up functions](#622-running-simulations-using-setup-files)
        - [6.2.3 Writing your matlab script from scratch](#623-writing-your-complete-matlab-script-yourself)
  - [7 MULTEM input parameters](#7-multem-input-parameters)
    - [7.1 Adding paths](#71-add-multem-paths-to-workspace)
    - [7.2 Setting system configuration](#72-specify-system-configuration)  
    - [7.3 Input parameters](#73-specify-input-parameters)
      - [7.3.1 Default parameters](#731-set-up-input-parameters)
      - [7.3.2 Simulation type](#732-set-the-simulation-type)
      - [7.3.3 Specimen interacion model](#733-set-the-electron-specimen-interaction-model)
      - [7.3.4 Load model file](#734-load-a-model-file-see-making-models-for-details)
      - [7.3.5 Output frequency](#735-set-slicing-and-when-to-output-results)
      - [7.3.6 Phonon interaction](#736-set-electron-phonon-interaction)
      - [7.3.7 Potential sampling (resolution and number of beams)](#737-set-potential-sampling)
      - [7.3.8 Acceleration voltage and beam tilt](#738-set-acceleration-voltage-and-beam-tilt)
      - [7.3.9 Illumination model](#739-set-the-illumination-model)
      - [7.3.10 Condeser lens parameters](#7310-set-condenser-lens-parmeters)
      - [7.3.11 Temporal incoherence](#7311-set-the-defocus-spread-function-temporal-incoherence)
      - [7.3.12 Spatial incoherence](#7312-set-the-source-spread-function-spatial-incoherence)
      - [7.3.13 Defocus reference](#7313-specify-the-zero-defocus-reference)
      - [7.3.14 STEM Detectors](#7314-define-stem-detectors)
      - [7.3.15 STEM Scan pattern](#7315-define-stem-scan-pattern)
    - [7.4 Setup functions](#74-setup-functions)
      - [7.4.1 STEM](#741-stem_setupm)
      - [7.4.2 HRTEM](#742-hrtem_setupm)
      - [7.4.3 CBED](#743-cbed_setupm)
      - [7.4.4 EWRS](#744-ewrs_setupm)
## 1 Introduction
mul2py is made to simplify analysis of MULTEM simulation results by providing tools to convert ".mat" and ".ecmat" files to HyperSpy signals. This is very much a developing package and there are many bugs and possible pitfalls, so it should be used carefully.

In addition, mul2py provides tools for making/converting modes, making images with the atomic model overlaid the results, and to export images that can be converted to movies (in e.g. ImageJ). Several other useful tools and examples are also provided, such as scanning CBED scripts, jupyter notebooks for creating models, and sample slurm job scripts.

[Up](#content)

## 2 Installation
mul2py can be found on github, and the easiest way to install it is by using git. Alternatively, the package may be downloaded manually and extracted instead of running Step 1. In addition, other environment tools than `virtualenv` may be used, such as `conda`, in which case some of the commands below must be changed. However, this is not covered at the moment.

[Up](#content)

### 2.1 Step 1: Clone mul2py
Navigate to where wou want to "install" mul2py and run
```bash
git clone https://github.com/TEM-Gemini-Centre/mul2py.git
```

[Up](#content)

#### 2.1.1 To update mul2py:
If there are updates to mul2py, this must be updated by navigating to where mul2py is installed and running
```bash
git init
git fetch origin
git checkout origin/branch
```

[Up](#content)

### 2.2 Step 2: Activate/create environment
If installing on the cluser, make sure that the correct Python version (3.7.x) is loaded along with its dependencies, i.e. you must first run
```bash
module load GCCcore/8.2.0 Python/3.7.2
```
Now you can create a new environment by navigating to where you want the environment to be installed by running
```bash
virtualenv mul2py-env
```
This environment can then be activated by running
```bash
source <path-to-mul2py-env>/bin/activate
```
[Up](#content)

### 2.3 Step 3: Add mul2py to environment
Add mul2py to the (active) environment by navigate to the root of mul2py and running
```bash
pip install --editable .
```
Note that the final `.` is very important.

mul2py should now be installed in editable mode and be visible when the environment is activated.

[Up](#content)

## 3 Convert MULTEM results
The main use of mul2py is to convert MULTEM results. Because MULTEM works in MATLAB, output are usually given in MATLABs `.mat` files. Because the data can be quite large, MATLABs 7.3 version of the `.mat` should also be used. These files are a variant of the HDF format, and can be loaded into Python by the `h5py` pacakge.

HDF data files are structured in `Groups`, similar to folders. These Groups may contain more Groups, or datasets. Datasets may be either actual data, or references to data which is stored somewhere/somehow in the file. This is important, because the data from MULTEM consists of a mix of structs, cells, and arrays, which are stored differently in the files. A field in a `struct` is saved as a Group (and sub-`struct`s are stored as sub-Groups). Arrays and numbers are stored as direct Datasets, while cells are stored as Referenced Datasets.

The point of mul2py is to make it easier to load this data into Python and construct Hyperspy signals for easier data inspection and analysis.

mul2py prefers to get siulation result files in the "ecmat" format (pronounced easy-mat) - a custom pseudo fileformat that is really just a normal MATLAB .mat file with data stored in specific fields. This helps to read and load the file in a convenient way:
```Python
from mul2py.buildtools import make_signal
signal = make_signal('my-simulation-results.ecmat')
```
As an alternative, mul2py can also create signals from regular .mat files, but this will still require data to be stored in specific fields. Most importantly, an output file must have the field `'output_multislice'` at its root, and an input file must have the field `'input_multislice'` at its root. Then, a signal can be made by
```Python
from mul2py.buildtools import load_input, load_output
input = load_input('my-simulation-input-parameters.mat') #This returns a dictionary with all the input parameters
signal, output = load_output('my-simulation-output.mat') #This returns an uncalibrated signal and the contents of the loaded file. This will eat memory quite fast (it duplicates the data a few times!) for large simulations.
```
[Up](#content)

### 3.1 Working concept
mul2py contains the subpackage `mul2py.io` that provides a tool to load the contents of a `hdf` file into an object with dynamic fields. These objects are used by the functions in `mul2py.buildtools.builders` to create hyperspy signals from the data files.

The next part is not important for the use of this package, as the user should not need to directly access the hdf reader and contents and instead use builders to create their signal. However, if the user stores his/hers MULTEM data in a different way, he/she must create his/hers own builder based on the `HDFReader` and `HDFContent` objects.

`mul2py.io.HDFReader` is an object that reads in the contents of a hdf file into `mul2py.HDFContent` objects. A `HDFContent` object contains more `HDFContent` objects in its `HDFContent._content` attribute (which may be accessed by the getter `HDFContent.content`), and the separate fields of this content is also stored as separate attributes in the object. 
For example, if a HDF file "my-data.mat" contains the groups "Level1.Level2.Data" and "Level1.Level2.Metadata", this file can be read by `my_data = mul2py.io.hdf.HDFReader(my-data.mat)`, and the contents of various groups may be accessed by `my_data.Level1.Level2.Data.content` or `my_data.Level1.Level2.Data()`.
Specifically, if a MULTEM `input_multislice` struct is stored as `input.mat`, the acceleration voltage of the simulation, and the STEM detectors can be accessed by  
```Python
from mul2py.io.hdf import HDFContent, HDFReader
data = HDFReader("input.mat")
E_0 = data.input_multislice.E_0()
det0_outer_ang = data.input_multislice.detector.cir.outer_ang.outer_ang0.resolved_reference()
```
The last line in the previous code block require some detailed explanation. The detectors in MULTEM are stored as cells:
```MATLAB
input_multislice.detector.cir(1).outer_ang = 40;
``` 
Now, when this is loaded by h5py, the content of `input_multislice.detector.cir` will be an array of references to other datasets (this is just how MATLAB stores cells in HDF files, and we have to deal with the problems it causes). When unpacking these references, `mul2py` goes through all the references and resolves them into separate variables with a suffix. In this case, the outer angles of detector 1 and 2 are stored as `outer_ang.outer_ang0` and `outer_ang.outer_ang1`, respectively. Because these are resolved references, the data in each of these fields is stored in the `resolved_reference` attribute, so it should be accessed by calling `outer_ang.outer_ang0.resolved_reference()`. The final paranthesis is the mul2py way of getting the value of a `HDFContent` object.

[Up](#content)

## 4 Plotting images and exporting data for movies
Sometimes, when plotting images, it is also useful to add the atoms in the model to the image. mul2py provides tools for this through the `mul2py.exporttools` subpackage.
```Python
from mul2py.buildtools import make_signal
from mul2py.exporttools.image import make_image
signal = make_signal('my_simulation_results.ecmat')

# Plot the simulation results of the current inav of `signal` and mark the atoms above this slice (the part where the beam has travelled):
fig, ax = make_image(signal)

# Plot the simulation results of a given inav of `signal`:
fig, ax = make_image(signal, inav=-1)

# Plot the simulation results of the current inav of `signal` and change the marker appearance of the atoms.
fig, ax = make_image(signal, markers={'cmap': 'jet'})
fig, ax = make_image(signal, markers={'edgecolors': 'k', 'facecolors': 'r', 's': 20, 'lw': 2, 'marker': 'x', 'alpha': 0.75})
```

It can also be very useful to export many images that can later be compiled to a video showing the scan pattern of the beam for instance:
```Python
import numpy as np
from mul2py.buildtools import make_signal
from mul2py.exporttools.movie import make_movie
signal = make_signal('mu_simulation_results.ecmat')
inavs = np.arange(0, len(signal)-1, 1) #Choose which slices to use in the movie. This line chooses all the slices - which is probably overkill in most cases (usually useful with certain thicknesses or more slices in the beginning and then fewer in the later parts of the stack.
make_movie('my_simulation_results', signal=signal, inavs = inavs)
```

[Up](#content)

## 5 Making models
Before MULTEM simulations can be run, you must make your atomistic model. MULTEM wants to get these models in a particular data format that includes e.g. the rms3D values (related to the Debye-Waller factors and temperature), and mul2py provides a way for converting .cif files to the matlab format that MULTEM uses:
```Python
from mul2py.io.models import save_multem_model
dwfs = {
    12: 1.8553,
    13: 0.7993,
    14: 0.4965,
    79: 0.6346,
} #The Debye-Waller factors of the atomic species in the model at a given temperature given as a dictionary {Z:B}, where Z is the atomic number and B is the debye-waller factor.
dz = 2.025 #The slice thickness in Å
save_multem_model('my-model.mat', 'my-model.cif', B = dwfs, dz = dz)
```

Alternatively, you can use [ASE](https://wiki.fysik.dtu.dk/ase/) to create your models directly. First, load a premade unit cell and duplicate it to a bigger slab, then make whatever changes you want, and then export it in a format that MULTEM likes by using mul2py.
```Python
import numpy as np

from mul2py.io.models import save_multem_model

from ase.visualize import view
from ase.io import read

#Make a pure Al slab
na, nb, nc = 10, 10, 20 #make a slab na * nb * nc unit cells large
Al = read('Al.cif') #Load a premade Al unit cell cif file.
a, b, c, alpha, beta, gamma = Al.get_cell_lengths_and_angles() #Get the lattice parameters of the model

slab = Al*[na, nb, nc] #Duplicate the unit cell to a larger slab

# Next, make some make some changes (skip this part if you only want the pure model). In this case, add an "L" shape of gold atoms in the middle of the model.
new_symbol = 'Au'
z_min = min([atom.position[-1] for atom in slab]) #Only make changes from this height
z_max = max([atom.position[-1] for atom in slab]) #Stop when this height is reached
    
positions = np.array([
    (0, 0),
    (0, 1),
    (0, 2),
    (0, 3),
    (1, 0),
    (2, 0),
]) #Make changes at these atomic columns (will be shifted to the middle later)

h = np.max(positions[:, 0]) + 1 #The height of the shape
w = np.max(positions[:,1]) + 1 #The width of the shape

offset = (na - int(h/2), nb - int(w/2)) #The x-y offset the shape to bring it into the middle

positions = positions + offset #Ofset the positions
positions = positions * a/2 #Scale the positions to the lattice and get them in Å.

#Loop through all atoms and change the atoms at the given positions
_ = [[atom.set('symbol', new_symbol) for atom in slab if (abs(atom.position[0] - x) < 1E-2 and abs(atom.position[1] - y) < 1E-2 and (z_min<= atom.position[2] <= z_max))] for x, y in positions]

#Inspect the model
view(slab)

# Finally, save the model to a mat file:
dwfs = {
    12: 1.8553,
    13: 0.7993,
    14: 0.4965,
    79: 0.6346,
}

save_multem_model('my-model.mat', slab, dz = a/2, B = dwfs)
```

[Up](#content)

## 6 Running MULTEM simulations
MULTEM simulations are run in MATLAB by calling
```MATLAB
output_multislice = IL_MULTEM(system_conf, input_multislice);
```
The `system_conf` and `input_multislice` are `struct` variables that are used to define the system configuration (CPU/GPU) and the input to the multislice simulation (specifying acceleration voltage, mode, specimen/atoms etc). Different fields in `input_multislice` are used depending on what kind of simulation experiment you are doing (HRTEM, STEM, CBED, etc), and in some cases, different fields have different units. Because most simulations of a particular type are very similar in terms of e.g. aberrations, defocus spread functions, and so on, mul2py provides some setup functions to help users focus on the parameters that are important for their particular case (the potential sampling, scanning, and convergence angles for instance). These setup functions take some required and optional parameters for certain fields, and use sensible values for the other input parameters that users don't need to consider as much. Please see However, when performing simulations to be used in publications and as a foundation for research, the user should always make their `input_multislice` themselves to be sure exactly what they are doing. 

In addition to the setup functions, mul2py also comes with a set of functions to wrap the setup, simulation, and results creation into one handy command. This is useful, because even if the input parameters themselves can be generated easily with the setup functions, performing the simulations by running a script that defines the `system_conf` and results construction is cumbersome - especially if you want to convert the results right after in the same shell script. Therefore, the `run_<simulation-type>_simulation()` functions are very useful for running your simulations without the need for a separate matlab script and to pass on variables directly from e.g. shell scripts instead.

[Up](#content)

### 6.1 Running simulations from shell
Running simulations "directly" from a shell script requires the use of a `run_<simulation-type>_simulation()` function, for example by calling `matlab -nodisplay -nodesktop -nosplash -r "run_HRTEM_simulation('<model-name>.mat'`. These functions require some information you must pass on in the shell script, such as the model name, and other optional parametes such as the sampling resolution. These values can also be passed on to the function using normal shell or bash syntax. For a simulation on the IDUN cluster, a shell script could like like this:
```shell script
#!bin/bash

#SBATCH --partition=GPUQ #Which partition to use
#SBATCH --time=00-00:15:00 #The allocated time
#SBATCH --job-name="Si2HfHRTEM" #The job name
#SBATCH --output=Si2Hf_HRTEM-%A.out #Where to put screen-output
#SBATCH --nodes=1 #The number of nodes
#SBATCH --ntasks-per-node=4 #The number of CPUs per node
#SBATCH --mem=64000 #The memory you need
#SBATCH --account=share-nv-fys-tem #Which account to use
#SBATCH --gres=gpu:1 #Whether or not you require GPU resources. Comment out if you are running on CPUQ

#Print some info about the job
echo "we are running from this directory: $SLURM_SUBMIT_DIR"
echo "The name of the job is: $SLURM_JOB_NAME"

echo "The job ID is $SLURM_JOB_ID"
echo "The job was run on these nodes: $SLURM_JOB_NODELIST"
echo "Number of nodes: $SLURM_JOB_NUM_NODES"
echo "We are using $SLURM_CPUS_ON_NODE cores"

echo "We are using $SLURM_CPUS_ON_NODE cores per node"
echo "Total of $SLURM_NTASKS cores"

#Load required modules
module load foss/2016a
module load CUDA/8.0.61
module load MATLAB/2017a

#*** Define some parameters to be passed to the run_<simulation-type>_simulation() function ***

#device=1 #CPU
device=2 #GPU

nx=1024 #The potential sampling to use (in X)
ny=1024 #The potential sampling to use (in Y)
instrument="2100F" #The instrument to use (determines aberrations, and will later determine source spread functions and defocus spread functions as well). Should be "2100F" or "ARM200F".

simulation_name="Si2Hf_101_8x3x150_HRTEM_${nx}x${ny}" #Name of your simulation. Results will be named "${simulation_name}_results.ecmat". Should match that of the model in this script
model_name="${simulation_name}.mat" #The model name
simulation_name="${simulation_name}_${instrument}" #Append the instrument to the simulation name - useful for performing simulations of the same model with different instruments

echo "Running HRTEM simulation: ${simulation_name}"

#Run a HRTEM simulation
# This line performs the MATLAB stuff. Remember to put `' '` around variables that should be passed as strings!
matlab -nodisplay -nodesktop -nosplash -r "addpath('/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py/mul2py/matlab'); run_HRTEM_simulation('$model_name','nx', $nx, 'ny', $ny, 'instrument', '$instrument', 'simulation_name', '${simulation_name}', 'print_parser', 1, 'print_details', 1, 'cpu_nthreads', $SLURM_NTASKS, 'device', $device);"

echo "Converting results to HyperSpy format using mul2py"
module load GCCcore/.8.2.0 Python/3.7.2
source /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py-env/bin/activate
python /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py/mul2py/examples/convert_ecmat.py "${simulation_name}_results.ecmat"

scontrol show job ${SLURM_JOB_ID} -d
```

[Up](#content)

### 6.2 Calling MATLAB script from shell
Alternatively, the simulation can be performed in a separate MATLAB script (called e.g. "multem.m") that is then run from the shell script with `matlab -nodisplay -nodesktop -nosplash -r "multem.m"`. In this case, there are many options for how to create your matlab script. You can either have a matlab script that does everything from defining input parameters, running the simulation, and constructing your results, or you can have a matlab script that uses some or all of the help functions of mul2py. This makes the shell script appear cleaner, but you must keep track of your input arguments and what your outputs are called yourself. A typical shell script would therefore look like this
```shell script
#!bin/bash

#SBATCH --partition=GPUQ #Which partition to use
#SBATCH --time=00-00:15:00 #The allocated time
#SBATCH --job-name="Si2HfHRTEM" #The job name
#SBATCH --output=Si2Hf_HRTEM-%A.out #Where to put screen-output
#SBATCH --nodes=1 #The number of nodes
#SBATCH --ntasks-per-node=1 #The number of CPUs per node. Make sure this matches the number of CPUs in the matlab script!
#SBATCH --mem=64000 #The memory you need
#SBATCH --account=share-nv-fys-tem #Which account to use
#SBATCH --gres=gpu:1 #Whether or not you require GPU resources. Comment out if you are running on CPUQ. Make sure this matches the device you choose in the matlab script (device 1 or device 2)!

#Print some info about the job
echo "we are running from this directory: $SLURM_SUBMIT_DIR"
echo "The name of the job is: $SLURM_JOB_NAME"

echo "The job ID is $SLURM_JOB_ID"
echo "The job was run on these nodes: $SLURM_JOB_NODELIST"
echo "Number of nodes: $SLURM_JOB_NUM_NODES"
echo "We are using $SLURM_CPUS_ON_NODE cores"

echo "We are using $SLURM_CPUS_ON_NODE cores per node"
echo "Total of $SLURM_NTASKS cores"

#Load required modules
module load foss/2016a
module load CUDA/8.0.61
module load MATLAB/2017a

echo "Running HRTEM simulation: ${simulation_name}"

#Run a HRTEM simulation
# This line performs the MATLAB stuff. Make sure that the resource allocation (the SBATCH stuff) matches what you use in the matlab script.
matlab -nodisplay -nodesktop -nosplash -r "multem.m"

echo "Converting results to HyperSpy format using mul2py"
module load GCCcore/.8.2.0 Python/3.7.2
source /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py-env/bin/activate
#Make sure the output from the matlab script is called "HRTEM_results.ecmat"
python /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py/mul2py/examples/convert_ecmat.py "HRTEM_results.ecmat"

scontrol show job ${SLURM_JOB_ID} -d
```
Note that you must keep track of what your simulation results are called yourself, and that the resources you ask for in the shell script matches those you use in the matlab script.

[Up](#content)

#### 6.2.1 Running `run_<simulation-type>_simulation()` through a script
A MATLAB script using the `run_<simulaiton-type>_simulation()` function would look like this:
```matlab
addpath('/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py/mul2py/matlab');

run_HRTEM_simulation('Si2Hf_101_8x3x150_HRTEM_1024_1024.mat','nx', 1024, 'ny', 1024, 'instrument', '2100F', 'simulation_name', 'HRTEM', 'print_parser', 1, 'print_details', 1, 'cpu_nthreads', 1, 'device', 2);
```

[Up](#content)

#### 6.2.2 Running simulations using setup files
You can write longer matlab scripts to execute to have more control yourself, but this will also increase the chance of errors. You can still use the setup files as a basis, and then modify them in the script. For instance, you can create your "multem.m" script using `HRTEM_setup()` to prepare the simulation parameters, and then change e.g. the defocus value. However, this also require you to set up other things as well, such as the system configuration, the paths used by the script, and to make sure the output is stored in a sensible way.
```matlab
%%
clear all
clc
%% System configuration

system_conf.precision = 1;               % eP_Float = 1, eP_double = 2
system_conf.device = 2;                  % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 5; 			 % Does the number of CPU threads matter when running on GPU?
system_conf.gpu_device = 1;				 % Which GPU device to use (if several are registered)

%% Timestamp
start_time = datetime('now','TimeZone','local');
fprintf("Starting simulation script at %s\n", start_time);

%% Paths
MULTEM_path = "/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM";  % Path to MULTEM installation on the cluster
addpath(char(sprintf("%s/crystalline_materials", MULTEM_path)));        % Add the crystalline materials function to the path (used for making models on the go if you want)
addpath(char(sprintf("%s/matlab_functions", MULTEM_path)));             % Add MULTEM matlab functions, such as `multem_default_values()` to the path.
addpath(char(sprintf("%s/mex_bin", MULTEM_path)));                      % Add the core MULTEM stuff to run simulations

addpath(char("/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py/mul2py/matlab")) % Add mul2py matlab scripts/functions to path

%% output_details
simulation_name = "HRTEM";
output_path = ".";
mkdir(char(output_path));

%% Load simulation parameters. `MULTEM_input.mat` should contain a struct called `input_multislice` with all relevant simulation parameters given in its fields, including the atomistic model.
input_multislice = HRTEM_setup("/lustre1/work/emilc/AlSiMgHf/MULTEM/101/Si2Hf_101_8x3x150_HRTEM_1024.mat", "nx", 1024, "ny", 1024, "bwl", 0, "phonons", 20, "instrument", "2100F", "multem_path", MULTEM_path);
input_multislice.obj_lens_c_10 = scherzer_defocus(input_multislice.E_0, input_multislice.obj_lens_c_30);

%% Run simulation
clear il_MULTEM;
tic;
output_multislice = il_MULTEM(system_conf, input_multislice);
toc;

%% Create results struct
results.input = input_multislice;
results.system = system_conf;
results.images = zeros(input_multislice.nx, input_multislice.ny, length(output_multislice.data));

%% Fill results structure
if length(output_multislice.data) == 1
    results.images(:,:, 1) = transpose(output_multislice.data.m2psi_tot);
else
    for t = 1:length(output_multislice.data)
        results.images(:, :, t) = transpose(output_multislice.data(t).m2psi_tot);
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

#### 6.2.3 Writing your complete matlab script yourself
Writing your own MATLAB script for doing MULTEM simulations is not difficult and is the normal way of doing it, but it requires some experience and trial-and-error. In the next sections, the steps are gone through step-by-step, and the various fields are explained. However, for more detailed explanations, please see the official MULTEM distribution pages.

## 7 MULTEM input parameters
### 7.1 Add MULTEM paths to workspace
Before you can run MULTEM, you must add the MULTEM paths to the workspace. This can be done by the following MATLAB code
```MATLAB
MULTEM_path = "/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM";  % Path to MULTEM installation on the cluster
addpath(char(sprintf("%s/crystalline_materials", MULTEM_path)));        % Add the crystalline materials function to the path (used for making models on the go if you want)
addpath(char(sprintf("%s/matlab_functions", MULTEM_path)));             % Add MULTEM matlab functions, such as `multem_default_values()` to the path.
addpath(char(sprintf("%s/mex_bin", MULTEM_path)));                      % Add the core MULTEM stuff to run simulations
```

[Up](#content)

### 7.2 Specify system configuration
To run simulations, you must also tell MULTEM what kind of system configuration you want it to use (number of CPUs or if you want to use GPU). This is done by creating a `struct` and defining certain fields with specific names (these fields will be extracted later by MULTEM; so it is important that they are named correctly!). The following defines a system configuration that will use 4 CPUs.
```MATLAB
system_conf.system_conf.precision = 1;              % eP_Float = 1, eP_double = 2
system_conf.device = 1;                             % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 4;                        % number of CPUs/threads
system_conf.gpu_device = 0;                         % which GPU device to use (only if system_conf.device=2). This enables different available GPUs to be selected.
```

[Up](#content)

### 7.3 Specify input parameters
Next, the input parameters to the multislice simulation must be defined. This is also done in a `struct` with specific fields. This is very sensitive to bugs, as you really have to specify the correct names. Normally, you would load the default MULTEM parameters using the `multem_default_values()` function, and then overwrite/change the fields you want/need. However, this means, that if you have typo, e.g. `input_multislice.smimulation_type=11` instead of the correct `input_multislice.simulation_type=11`, MULTEM will use whatever value that is assigned to `input_multislice.simulation_type` in `multem_default_values()`. In other words, be careful when setting up your simulation parameters! mul2py also provides a set of help-functions for setting up simulations with common input parameters (and with easy modification of these), see section [6.4](#64-setup-functions) for a list of these functions and how to use them.

The following code sets up a STEM simulation in MULTEM. For other simulation types, other fields in `input_Multislice` are used and must be changed instead. The best way of figuring out what parameters to use/change for the different simulation modes is by chekcing out the MULTEM examples. 

[Up](#content)

#### 7.3.1 Set up input parameters
```MATLAB
input_multislice = multem_default_values();         % Loads default values
```

[Up](#content)

#### 7.3.2 Set the simulation type
```MATLAB
input_multislice.simulation_type = 11; % STEM type
```

[Up](#content)

#### 7.3.3 Set the electron-specimen interaction model
```MATLAB
input_multislice.interaction_model = 1;              % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_type = 6;                 % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6
```

[Up](#content)

#### 7.3.4 Load a model file? See "Making models" for details
```MATLAB
load(p.Results.model_path);                 % This should load `spec_atoms`, `spec_lx`, `spec_ly`, `spec_lz`, `spec_dz`, `a`, `b`, `c`, `na`, `nb`, and `nc` into the MATLAB workspace.
input_multislice.spec_atoms = spec_atoms;   % The atomic species, positions, and rms3d values (Debye-Waller factors). See "making models" for more details on how to make this
input_multislice.spec_lx = spec_lx;         % The model X-dimension
input_multislice.spec_ly = spec_ly;         % The model Y-dimension
input_multislice.spec_lz = spec_lz;         % The model Z-dimension
input_multislice.spec_dz = spec_dz;         % The slicing thickness
input_multislice.spec_cryst_a = a;          % The lattice a-parameter (not used directly)
input_multislice.spec_cryst_b = b;          % The lattice b-parameter (not used directly)
input_multislice.spec_cryst_c = c;          % The lattice c-parameter (not used directly)
input_multislice.spec_cryst_na = na;        % The number of unit cells in a (not used directly)
input_multislice.spec_cryst_nb = nb;        % The number of unit cells in b (not used directly)
input_multislice.spec_cryst_nc = nc;        % The number of unit cells in c (not used directly)
```

[Up](#content)

#### 7.3.5 Set slicing and when to output results
```MATLAB
input_multislice.potential_slicing = 2; %Slice the potential by slice projection
input_multislice.thick_type = 2; % Get output by slices
input_multislice.thick = (0:input_multislice.spec_dz:input_multislice.spec_lz-input_multislice.spec_dz); % Get outputs for these thicknesses
```

[Up](#content)

#### 7.3.6 Set electron-phonon interaction
```MATLAB
input_multislice.pn_model = 3;                       % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.pn_coh_contrib = 0;
input_multislice.pn_single_conf = 0;                 % 1: true, 0:false (extract single configuration)
input_multislice.pn_nconf = 20;                      % true: specific phonon configuration, false: number of frozen phonon configurations
input_multislice.pn_dim = 110;                       % phonon dimensions (xyz)
input_multislice.pn_seed = 300183;                   % Random seed(frozen phonon)
```

[Up](#content)

#### 7.3.7 Set potential sampling.
This affects resolution/maximum scattering angle and affects the simulation time alot.
```MATLAB
input_multislice.nx = 1024;     % number of pixels in x
input_multislice.ny = 1024;     % number of pixels in y
input_multislice.bwl = 1;   % Band-width limit, 1: true, 0:false
```

[Up](#content)

#### 7.3.8 Set acceleration voltage and beam tilt
```MATLAB
input_multislice.E_0 = 200.00;  % Acceleration Voltage (keV)
input_multislice.theta = 0.0;   % Tilt ilumination (deg)
input_multislice.phi = 0.0;     % Tilt ilumination (deg)
```

[Up](#content)

#### 7.3.9 Set the illumination model
```MATLAB
input_multislice.illumination_model = 1;        % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
input_multislice.temporal_spatial_incoh = 1;    % 1: Temporal and Spatial, 2: Temporal, 3: Spatial
```

[Up](#content)

#### 7.3.10 Set condenser lens parmeters
More aberrations can be defined, but C10 and C12 are the most important
```MATLAB
input_multislice.cond_lens_m = 0;                   % Vortex momentum
input_multislice.cond_lens_inner_aper_ang = 0.0;    % Inner aperture (mrad) 
input_multislice.cond_lens_outer_aper_ang = 27.42;  % Outer aperture (mrad) defines convergence semi-angle

input_multislice.obj_lens_c_10 = 0.00;  % Defocus [Å]
input_multislice.cond_lens_c12 = 0.00;  %Condenser lens spherical aberration [mm]
```

[Up](#content)

#### 7.3.11 Set the defocus spread function (temporal incoherence?)
Only used if the illumination model is set to 4 (numerical integration) I think.
```MATLAB
dsf_sigma = il_iehwgd_2_sigma(32);                  % from defocus spread to standard deviation
input_multislice.cond_lens_dsf_sigma = dsf_sigma;   % standard deviation (�)
input_multislice.cond_lens_dsf_npoints = 5;         % # of integration points. It will be only used if illumination_model=4
```

[Up](#content)

#### 7.3.12 Set the source spread function (spatial incoherence)
Only used if the illumination model is set to 4 (numerical integration) I think.
##### 7.3.12.1 Converged mode:
```MATLAB
ssf_sigma = il_hwhm_2_sigma(0.45); % half width at half maximum to standard deviation
input_multislice.cond_lens_ssf_sigma = ssf_sigma;  	% standard deviation: For parallel ilumination(�^-1); otherwise (�)
input_multislice.cond_lens_ssf_npoints = 4;         % # of integration points. It will be only used if illumination_model=4
```

[Up](#content)

##### 7.3.12.1 Parallel mode:
```MATLAB
ssf_sigma = il_mrad_2_sigma(input_multislice.E_0, 0.02); % mrad to standard deviation
input_multislice.cond_lens_ssf_sigma = ssf_sigma;  	% standard deviation: For parallel ilumination(�^-1); otherwise (�)
input_multislice.cond_lens_ssf_npoints = 4;         % # of integration points. It will be only used if illumination_model=4
```

[Up](#content)

#### 7.3.13 Specify the zero defocus reference
```MATLAB
input_multislice.cond_lens_zero_defocus_type = 1;   % eZDT_First = 1, eZDT_User_Define = 2
input_multislice.cond_lens_zero_defocus_plane = 0;
```

[Up](#content)

#### 7.3.14 Define STEM detectors
```MATLAB
input_multislice.detector.type = 1;  % eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3

% Set collection angles (in mrad) of circular detectors. You can add as many detectors as you want by defining `input_multislice.detector.cir(x).inner_ang = <some-angle>; input_multislice.detector.cir(x).outer_ang = <some-ang>`
input_multislice.detector.cir(1).inner_ang = 5.00;      % inner collection semi-ngle of detector 1 [mrad]
input_multislice.detector.cir(1).outer_ang = 40.00;     % outer collection semi-angle of detector 1 [mrad]
input_multislice.detector.cir(2).inner_ang = 48.00;     % inner collection semi-angle of detector 2 [mrad]
input_multislice.detector.cir(2).outer_ang = 203.00;    % outer collection semi-angle of detector 2 [mrad]
```

[Up](#content)

#### 7.3.15 Define STEM scan pattern
```MATLAB
input_multislice.scanning_ns = 25; % 25 probe positions
input_multislice.scanning_x0 = 0.00; % scan start in x [Å]
input_multislice.scanning_xe = 4.05; % scan stop in x [Å]
input_multislice.scanning_y0 = 0.00; % scan start in y [Å]
input_multislice.scanning_ye = 4.05; % scan stop in y [Å]
input_multislice.scanning_periodic = 0; % omit last scan row/column to enable periodic replication of images? 0: false, 1: true
```

[Up](#content)

### 7.4 Setup functions
As you can see, setting up a simulation requires quite many lines of code, and are very sensitive to typos and errors. Therefore, `mul2py` also provides some help functions to set up some usual simulation types:
 - `STEM_setup("model_path", convergence_angle, detectors)`
 - `HRTEM_setup("model_path")`
 - `CBED_setup("model_path", convergence_angle)`
 - `EWRS_setup("model_path", convergence_angle)`
 
 Each of these functions take several optional inputs to ease configuration of simulation files.

[Up](#content)

#### 7.4.1 STEM_setup.m
`STEM_setup.m` takes three required parameters and a number of optional parameters.

Example
 ```MATLAB
addpath(char(sprintf("%s/mul2py/mul2py/matlab/", path-to-mul2py))) %Make MATLAB aware of mul2py matlab functions

model_path = "path/to/my/model.mat";    %Path to my model

alpha = 27.42;                          %Convergence semi-angle in [mrad]

clear collection_angles
collection_angles(1).inner_ang = 0;     %Inner collection semi-angle of detector 1 (mrad)
collection_angles(1).outer_ang = 40;    %Outer collection semi-angle of detector 1 (mrad)
collection_angles(2).inner_ang = 48;    %Inner collection semi-angle of detector 2 (mrad)
collection_angles(2).outer_ang = 200;   %Outer collection semi-angle of detector 2 (mrad

input_multislice = STEM_setup(model_path, alpha, collection_angles)
```
 
 ##### Required parameters:
  - `model_path`: `str` or `char`
    - Path to your model
  - `alpha`: `float`
    - Convergence semi-angle in mrad
  - `collection_angles`: `struct`
    - Collection semi-angles in mrad
    
##### Optional parameters:
  - `defocus=0.00`: `float`
    - Defocus of condenser lens [Å]
  - `nx=1024`: `int`
    - Potential pixels in x
  - `ny=1024: `int`
    - Potential pixels in y
  - `bwl=1`: `int` or `float` or `bool`
    - Bandwidth limit?
  - `center_x=nan`: `float`
    - Center of scan pattern along x [Å]. If `nan`, defaults to `spec_lx/2`
  - `center_y=nan`: `float`
    - Center of scan pattern along y [Å]. If `nan`, defaults to `spec_ly/2`
  - `scanning_ns=10`: `int`
    - Number of scanning points
  - `scan_width=nan`: `float`
    - Width of scan (x dimension) [Å]. If `nan`, defaults to `spec_cryst_a`
  - `scan_height=nan`: `float`
    - Height of scan (y dimension) [Å]. If `nan`, defaults to `spec_cryst_b`
  - `E0=200.00`: `float`
    - Acceleration voltage [kV]
  - `phonons=20`: `int`
    - Number of phonon configurations
  - `thick_type=2`: `int`
    - When to output results (through the thickness or after each slice)
  - `thicknesses=0`: `int` or `array`
    - Output at these thicknesses if `thick_type=2`. If `0`, defaults to each slice.
  - `instrument=""`: `str` or `char`
    - Which instrument? Decides the aberrations. If `""` sets all aberrations to 0. Other valid instruments are `"ARM200F"` and `"2100F"` which will set aberrations to some sensible values for these instruments
  - `print_parser=0`: `int`, `float`, or `bool`
    - Print the parameters passed to this function?
  - `print_details=1`: `int`, `float`, or `bool`
    - Print details of the simulation (wavelength, resolution, maximum scattering angle etc)
  - `MULTEM_path="/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM"`: `str` or `char`
    - The path to MULTEM. Defaults to IDUN project folder.

[Up](#content)

#### 7.4.2 HRTEM_setup.m
`HRTEM_setup.m` takes one required parameter and a number of optional parameters.

Example
 ```MATLAB
addpath(char(sprintf("%s/mul2py/mul2py/matlab/", path-to-mul2py))) %Make MATLAB aware of mul2py matlab functions

model_path = "path/to/my/model.mat";    %Path to my model

input_multislice = HRTEM_setup(model_path)
```
 
 ##### Required parameters:
  - `model_path`: `str` or `char`
    - Path to your model
        
##### Optional parameters:
  - `mode="converged"`: `str` or `char`
    - Mode of beam. Should be either `"parallel"` or `"converged"`
  - `defocus=0.00`: `float`
    - Defocus of objective lens [Å]
  - `nx=1024`: `int`
    - Potential pixels in x
  - `ny=1024: `int`
    - Potential pixels in y
  - `bwl=1`: `int` or `float` or `bool`
    - Bandwidth limit?
  - `E0=200.00`: `float`
    - Acceleration voltage [kV]
  - `phonons=20`: `int`
    - Number of phonon configurations
  - `thick_type=2`: `int`
    - When to output results (through the thickness or after each slice)
  - `thicknesses=0`: `int` or `array`
    - Output at these thicknesses if `thick_type=2`. If `0`, defaults to each slice.
  - `instrument=""`: `str` or `char`
    - Which instrument? Decides the aberrations. If `""` sets all aberrations to 0. Other valid instruments are `"ARM200F"` and `"2100F"` which will set aberrations to some sensible values for these instruments
  - `print_parser=0`: `int`, `float`, or `bool`
    - Print the parameters passed to this function?
  - `print_details=1`: `int`, `float`, or `bool`
    - Print details of the simulation (wavelength, resolution, maximum scattering angle etc)
  - `MULTEM_path="/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM"`: `str` or `char`
    - The path to MULTEM. Defaults to IDUN project folder.

[Up](#content)

#### 7.4.3 CBED_setup.m
`CBED_setup.m` takes two required parameters and a number of optional parameters.

Example
 ```MATLAB
addpath(char(sprintf("%s/mul2py/mul2py/matlab/", path-to-mul2py))) %Make MATLAB aware of mul2py matlab functions

model_path = "path/to/my/model.mat";    %Path to my model

input_multislice = CBED_setup(model_path)
```
 
 ##### Required parameters:
  - `model_path`: `str` or `char`
    - Path to your model
  - `alpha`: `float`
    - Convergence semi-angle of beam.
        
##### Optional parameters:
  - `defocus=0.00`: `float`
    - Defocus of condenser lens [Å]
  - `nx=1024`: `int`
    - Potential pixels in x
  - `ny=1024: `int`
    - Potential pixels in y
  - `bwl=1`: `int` or `float` or `bool`
    - Bandwidth limit?
  - `E0=200.00`: `float`
    - Acceleration voltage [kV]
  - `phonons=20`: `int`
    - Number of phonon configurations
  - `thick_type=2`: `int`
    - When to output results (through the thickness or after each slice)
  - `thicknesses=0`: `int` or `array`
    - Output at these thicknesses if `thick_type=2`. If `0`, defaults to each slice.
  - `instrument=""`: `str` or `char`
    - Which instrument? Decides the aberrations. If `""` sets all aberrations to 0. Other valid instruments are `"ARM200F"` and `"2100F"` which will set aberrations to some sensible values for these instruments
  - `print_parser=0`: `int`, `float`, or `bool`
    - Print the parameters passed to this function?
  - `print_details=1`: `int`, `float`, or `bool`
    - Print details of the simulation (wavelength, resolution, maximum scattering angle etc)
  - `MULTEM_path="/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM"`: `str` or `char`
    - The path to MULTEM. Defaults to IDUN project folder.

[Up](#content)

#### 7.4.4 EWRS_setup.m
`EWRS_setup.m` takes two required parameters and a number of optional parameters.

Example
 ```MATLAB
addpath(char(sprintf("%s/mul2py/mul2py/matlab/", path-to-mul2py))) %Make MATLAB aware of mul2py matlab functions

model_path = "path/to/my/model.mat";    %Path to my model

input_multislice = EWRS_setup(model_path)
```
 
 ##### Required parameters:
  - `model_path`: `str` or `char`
    - Path to your model
  - `alpha`: `float`
    - Convergence semi-angle of beam.
        
##### Optional parameters:
  - `mode="converged"`: `str` or `char`
    - Mode of beam. Should be either `"parallel"` or `"converged"`
  - `defocus=0.00`: `float`
    - Defocus of condenser lens [Å]
  - `nx=1024`: `int`
    - Potential pixels in x
  - `ny=1024: `int`
    - Potential pixels in y
  - `bwl=1`: `int` or `float` or `bool`
    - Bandwidth limit?
  - `E0=200.00`: `float`
    - Acceleration voltage [kV]
  - `phonons=20`: `int`
    - Number of phonon configurations
  - `thick_type=2`: `int`
    - When to output results (through the thickness or after each slice)
  - `thicknesses=0`: `int` or `array`
    - Output at these thicknesses if `thick_type=2`. If `0`, defaults to each slice.
  - `instrument=""`: `str` or `char`
    - Which instrument? Decides the aberrations. If `""` sets all aberrations to 0. Other valid instruments are `"ARM200F"` and `"2100F"` which will set aberrations to some sensible values for these instruments
  - `print_parser=0`: `int`, `float`, or `bool`
    - Print the parameters passed to this function?
  - `print_details=1`: `int`, `float`, or `bool`
    - Print details of the simulation (wavelength, resolution, maximum scattering angle etc)
  - `MULTEM_path="/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM"`: `str` or `char`
    - The path to MULTEM. Defaults to IDUN project folder.
    
[Up](#content)
