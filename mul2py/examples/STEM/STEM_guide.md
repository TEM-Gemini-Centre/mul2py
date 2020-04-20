# STEM guide
This is a guide on how to prepare, setup and run a STEM MULTEM simulation on the IDUN cluster on NTNU

## General prepartions
The first things you must consider when performing any kind of STEM simulation is
  - Which high-tension do you want to use?
  - What convergence angle do you want to use?
  - Which collection angles do you need?
    - Doesn't affect simulation time that much - but can influence finished file size
  - Aberrations and defocus?
  - How many phonon configurations do you need?
    - Usually safe to use about 20
  - Simulation resolution (NOT scan-resolution)?
    - Determined by the potential sampling (number of beams)
    - Large impact on simulation time
  - How many scan pixels do you need?
    - Large impact on simulation time
  
The answer to many of these questions will influence how you make your model and how you should run your simulation. For instance, if you use too many potential pixels, phonon configurations, and scan pixels, your simulation will be too timeconsuming, even for IDUN. In that case, you can divide your simulation into several simulations with fewer scan pixels and stich them together afterwards, if you really need that large simulations.

## Making your model
After deciding on these questions, you must make your model. You should take a look at the example notebook in [https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/examples/Models](https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/examples/Models) to get a more thorough introduction to making models for MULTEM. The important thing here, is that your model size will decide (together with the potential sampling) on the maximum scattering angle you will get - and this must match your collection angles for your detectors!

## Preparing your simulation
To make the MATLAB script for running MULTEM simulations, it is probably easiest to use `mul2py`s STEM setup function: [https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/matlab/STEM_setup.m](https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/matlab/STEM_setup.m). This function will set up most of what you need, and make it easier to get going. However, if you plan on using these simulations in a publication, you should prepare your simulation yourself (based on the example scripts provided by MULTEM).

### Step 1 - system configuration
You must tell MULTEM on what system resources you want to use. This is done by defining the following struct:
```MATLAB
system_conf.precision = 1;                           % eP_Float = 1, eP_double = 2
system_conf.device = 2;                              % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 1;
system_conf.gpu_device = 1;
```
Usually, floating point precision is fine (`system_conf.precision = 1`). Simulations run faster when using GPUs on IDUN (`system_conf.device=2`), in which case you will only need one CPU (`system_conf.cpu_nthread = 1`). Finally, you must choose which GPU you want to use (on IDUN there are two identical GPUs on each GPU node, so just set `system_conf.gpu_devide=1`).

### Step 2 - set up paths
It is often a good idea to set up paths at the beginning of your script (and maybe also try to write to the paths so that  you know they exist and that you have writing permissions):
```MATLAB
MULTEM_path = "/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM";  % Path to MULTEM installation on the cluster
addpath(char(sprintf("%s/crystalline_materials", MULTEM_path)));        % Add the crystalline materials function to the path (used for making models on the go if you want)
addpath(char(sprintf("%s/matlab_functions", MULTEM_path)));             % Add MULTEM matlab functions, such as `multem_default_values()` to the path.
addpath(char(sprintf("%s/mex_bin", MULTEM_path)));                      % Add the core MULTEM stuff to run simulations

%% output_details
simulation_name = "STEM";
output_path = "/lustre1/work/emilc/MULTEM/Test/"; %Path to put results
mkdir(output_path);
```
The first three `addpath()` statements are necessary for MATLAB to find MULTEM, and should always be included in your MULTEM scripts. The `MULTEM_path` variable in this example points to the path of the MULTEM installation on the TEM Gemini project folder on IDUN, but should be changed to your own path if you are running locally. The final three lines prepares the folder for your output and defines a name for your output files. You could also add some lines where you create and write a log file or something to test that you are allowed to write to the output folder.

### Step 3 - set up your simulation parameters
Now you can start making your `input_multislice` struct that defines the parameters used in the simulation. In this guide, we are using the `STEM_setup.m` function of `mul2py` to do this (you can change everything afterwards if you like):
```MATLAB
clear collection_angles
collection_angles(1).inner_ang = 0; %semi-angle (mrad)
collection_angles(1).outer_ang = 40; %semi-angle (mrad)
collection_angles(2).inner_ang = 48; %semi-angle (mrad)
collection_angles(2).outer_ang = 200; %semi-angle (mrad)

convergence_angle = 27; %semi-angle (mrad)

input_multislice = STEM_setup("Al_10x10x20.mat", convergence_angle, collection_angles, "phonons", 20, "nx", 1024, "ny", 1024, "instrument", "ARM200F", "multem_path", multem_path);
```
First, we define the collection angles we want to use as a cell-array with struct fields (inner and outer angles). Then, we define our convergence angle, before we call the setup function. This function takes three required parameters: The model path, the convergence angle, and the collection angles. Additionally, you can provide several optional parameters, such as the number of phonons, the potential sampling, the instrument (predefined aberration values), and the path to the MULTEM installation, among many more. You can also specify the scan start and stop values, but this is best done afterwards.

### Step 4 - Define your scan area
After defining your `input_multislice` struct, you can make changes to its fields (and add new ones). For instance, you can change the scan area of your simulation:
```MATLAB
input_multislice.scanning_ns = 25;
input_multislice.scanning_periodic = 0;
input_multislice.scanning_x0 = input_multislice.spec_lx/2 - input_multislice.spec_cryst_a/2;
input_multislice.scanning_xe = input_multislice.spec_lx/2 + input_multislice.spec_cryst_a/2;

input_multislice.scanning_y0 = input_multislice.spec_ly/2 - input_multislice.spec_cryst_b/2;
input_multislice.scanning_ye = input_multislice.spec_ly/2 + input_multislice.spec_cryst_b/2;
```
This will scan 25 pixels (`input_multislice.scanning_ns=25`) and include the last row and column (`input_multislice.scanning_periodic=0`). It will start half a unit cell from the center of the model (`input_multislice.scanning_x0 = input_multislice.spec_lx/2 - input_multislice.spec_cryst_a/2`) and scan one unit cell to end up at half a unit cell away from the centre of the model in the other direction (`input_multislice.scanning_x0 = input_multislice.spec_lx/2 + input_multislice.spec_cryst_a/2`).

### Step 5 - Execute simulation
When all is prepared, you can run the simulation by sending `system_conf` and `input_multislice` to the `il_MULTEM` function:
```MULTEM
clear il_MULTEM;
tic;
output_multislice = il_MULTEM(system_conf, input_multislice); 
toc;
```

### Step 6 - Store output
Output can be stored simply by `save("output_multislice.mat", "output_multislice", "-v7.3")`, or you can store them in a new struct to make it possible for `mul2py` to generate Hyperspy signals:
```MATLAB
%% Store input parameters in results struct
results.input = input_multislice;
results.system = system_conf;

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

%% Save the data
save(sprintf("%s/%s_results.ecmat", output_path, simulation_name), "results", "-v7.3");
```

## Running your simulation
When you have prepared your files (the model and the matlab script), upload them to the cluster, write a suitable SLURM script, and queue the job. It is important that the SLURM script is tuned to your simulation (run time, number of cores, GPU assignment, etc) and you should make this yourself. However, the following job script should work fine for this example and serve as a template. It requires to be in the same folder as the MULTEM MATLAB script called "STEM.m" in this case.

```shell script
#!/bin/bash

#SBATCH --partition=GPUQ
#SBATCH --time=01-20:0:00
#SBATCH --job-name="STEM"
#SBATCH --output=STEM-%A.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -- gres=gpu:1
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

module load foss/2016a
module load CUDA/8.0.61
module load MATLAB/2017a

echo "Running STEM simulation"
matlab -nodisplay -nodesktop -nosplash -r "STEM"

scontrol show job ${SLURM_JOB_ID} -d
```

## Converting your results
You can also use the cluster to convert your results through the `mul2py` package. The project folder on IDUN should include a virtual environment where `mul2py` is installed, and you can run the python script `convert_ecmat.py` to make a HyperSpy file from the output of the simulation.

```shell script
#!/bin/bash

#SBATCH --partition=CPUQ
#SBATCH --time=00-01:0:00
#SBATCH --job-name="STEMConvert"
#SBATCH --output=STEMConvert-%A.out
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

echo "Converting results to HyperSpy format using mul2py"
module load GCCcore/.8.2.0 Python/3.7.2
source /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py-env/bin/activate
python /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py/mul2py/examples/convert_ecmat.py /lustre1/work/emilc/MULTEM/Test/STEM_results.ecmat

scontrol show job ${SLURM_JOB_ID} -d
```
In this case, the shell script requires that the output from the simulation is stored at "/lustre1/work/emilc/MULTEM/Test/STEM_results.ecmat", and that the path to the python script "convert_ecmat.py" is correct. 

### About `convert_ecmat.py`
This script takes one argument, namely the path to a ".ecmat" file (output from the MULTEM simulation), creates a complete and scaled HyperSpy signal with sensible metadata and saves this signal to a file. The filename of the output file is the same as the input file and the extension is ".hspy". The script also takes one positional argument "--overwrite" that determines if any existing file will be overwritten (default is `True`). In other words, if you call `python convert_ecmat.py STEM_results.ecmat`, a file with the filename "STEM_results.hspy" will be created or overwritten. If you call `python convert_ecmat.py STEM_results.ecmat --overwrite False`, the same file will be created unless a file named "STEM_results.hspy" already exists in the directory.

## Postprocessing
After the output has been created (and possibly converted), it is time to make sense of it. If simply the `output_multislice` struct is written to file, it is best to do the postprocessing in MATLAB, or to export the data to some other format to treat them in your favourite tool. If an ".ecmat" file was created instead and has been converted to a ".hspy" file, you can use HyperSpy to postprocess your data.

### HyperSpy
To analyse your hyperspy STEM simulation signal, you should first load the data, and then do whatever analysis you require. One possible and useful method is to create series of virtual dark field images of particular atomic columns and see how the scattering from this column develops through the thickness. If your model has an atomic column situated at $x=0$ Å and $y=0$ Å, tou can make a thickness profile of that column by integrating the scattering intensity from its center out to about $1/4$ lattice parameter ($1.025$ Å) as in the following example:

```Python
import hyperspy.api as hs

signal = hs.load("STEM_results.hspy") #Load the converted ".ecmat" file
print(signal.metadata) #Print the metadata of the signal
print(signal.axes_manager) #Print some info about the axes

signal.plot() #Plot the signal
roi = hs.roi.CircleROI(cx=0, cy=0, r=1.025) #Make a region of interest. Make sure that the position of the roi is within the x-y-space of the signal
roi.add_widget(signal, axes=['x', 'y']) #Connect the roi to the signal
cropped_signal = roi(signal) #Extract the roi from the signal (does not affect the original signal)

scattering_thickness_profile = sum(cropped_signal, axes=['x', 'y']) #Sum the signal in the x and y axes

detector = 1 #Choose a detector in the signal
scattering_thickness_profile.inav[:, detector].plot() #Plot the thickness profile of the image stack for the selected detector
```  
For a more detailed explanation and guide, please see the notebook at [https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/examples/STEM/STEM_postprocessing.ipynb](https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/examples/STEM/STEM_postprocessing.ipynb).

## Complete example

The files you need to run a complete STEM simulation with MULTEM on the IDUN cluster is already provided in this folder. You can simply copy these files ("Al_10x10x20.mat", "STEM.m", and "STEM.slurm") to somewhere on your IDUN directory(e.g. "/lustre1/work/emilc/") and run the following command on IDUN to perform the simulation:

```bash
cd /lustre1/work/emilc/
sbatch STEM.slurm
```

To understand these files, continue reading

### The model file
The model file was made by the following Python code:
```Python
from ase.io import read
import mul2py as m2p

Al = read('Al.cif') #Load the crystal information file

na, nb, nc = 10, 10, 20 # Number of unit cells along a, b, and c crystal axes

slab = Al*[na, nb, nc] #Duplicate the model to make a slab

slab.center(axis=(0, 1)) #Center the slab in the x-y plane - leave it unchanged in z-direction

dwfs = {13: 0.1006}
m2p.io.save_multem_model('Al_10x10x20.mat', slab, B=dwfs)
```

### The simulation script
The somulation script "STEM.m" contains the following (or similar) MATLAB code
```MATLAB
%%
clear all
clc

%% System configuration
system_conf.precision = 1;                           % eP_Float = 1, eP_double = 2
system_conf.device = 2;                              % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 1; 			 % Does the number of CPU threads matter when running on GPU? EXPERIMENT!!
system_conf.gpu_device = 1;				 % MULTEM can only use one GPU device at the time? Only ask for a single GPU from IDUN, and use this.

%% Timestamp
start_time = datetime('now','TimeZone','local');
fprintf("Starting simulation script at %s\n", start_time);

%% Paths
MULTEM_path = "/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM";  % Path to MULTEM installation on the cluster
addpath(char(sprintf("%s/crystalline_materials", MULTEM_path)));        % Add the crystalline materials function to the path (used for making models on the go if you want)
addpath(char(sprintf("%s/matlab_functions", MULTEM_path)));             % Add MULTEM matlab functions, such as `multem_default_values()` to the path.
addpath(char(sprintf("%s/mex_bin", MULTEM_path)));                      % Add the core MULTEM stuff to run simulations

%% output_details
simulation_name = "STEM";
output_path = "/lustre1/work/emilc/MULTEM/Test/"; %Path to put results
mkdir(output_path);

%% Load simulation parameters. `MULTEM_input.mat` should contain a struct called `input_multislice` with all relevant simulation parameters given in its fields, including the atomistic model.
clear collection_angles
collection_angles(1).inner_ang = 0; %semi-angle (mrad)
collection_angles(1).outer_ang = 40; %semi-angle (mrad)
collection_angles(2).inner_ang = 48; %semi-angle (mrad)
collection_angles(2).outer_ang = 200; %semi-angle (mrad)

convergence_angle = 27; %semi-angle (mrad)

input_multislice = STEM_setup("Al_10x10x20.mat", convergence_angle, collection_angles, "phonons", 20, "nx", 1024, "ny", 1024, "instrument", "ARM200F", "multem_path", MULTEM_path);

%% Adjust scan pattern
input_multislice.scanning_ns = 25;
input_multislice.scanning_periodic = 0;
input_multislice.scanning_x0 = input_multislice.spec_lx/2 - input_multislice.spec_cryst_a/2;
input_multislice.scanning_xe = input_multislice.spec_lx/2 + input_multislice.spec_cryst_a/2;

input_multislice.scanning_y0 = input_multislice.spec_ly/2 - input_multislice.spec_cryst_b/2;
input_multislice.scanning_ye = input_multislice.spec_ly/2 + input_multislice.spec_cryst_b/2;

%% Run simulation
clear il_MULTEM;
tic;
output_multislice = il_MULTEM(system_conf, input_multislice); 
toc;

%% Store input parameters in results struct
results.input = input_multislice;
results.system = system_conf;

%% Construct results images
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

%% Save the data
save(sprintf("%s/%s_results.ecmat", output_path, simulation_name), "results", "-v7.3");
```

### The SLURM job script
And finally, the slurm job shell script "STEM.slurm" contains the following bash code
```shell script
#!/bin/bash

#SBATCH --partition=GPUQ
#SBATCH --time=01-20:0:00
#SBATCH --job-name="STEM"
#SBATCH --output=STEM-%A.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64000
#SBATCH --gres=gpu:1
#SBATCH --account=share-nv-fys-tem

echo "we are running from this directory: $SLURM_SUBMIT_DIR"
echo "The name of the job is: $SLURM_JOB_NAME"

echo "The job ID is $SLURM_JOB_ID"
echo "The job was run on these nodes: $SLURM_JOB_NODELIST"
echo "Number of nodes: $SLURM_JOB_NUM_NODES"
echo "We are using $SLURM_CPUS_ON_NODE cores"

echo "We are using $SLURM_CPUS_ON_NODE cores per node"
echo "Total of $SLURM_NTASKS cores"

module load foss/2016a
module load CUDA/8.0.61
module load MATLAB/2017a

echo "Running STEM simulation"
matlab -nodisplay -nodesktop -nosplash -r "STEM"

echo "Converting results to HyperSpy format using mul2py"
module load GCCcore/.8.2.0 Python/3.7.2
source /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py-env/bin/activate
python /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py/mul2py/examples/convert_ecmat.py STEM_results.ecmat

scontrol show job ${SLURM_JOB_ID} -d
```

## Postprocessing
Postprocess your data as you would with other TEM data, e.g. by making virtual images or thickness profile plots, e.g in a python script or a jupyter notebook.

### Spatial incoherence
If you have not run the simulation with spatial incoherence and numerical integration options (i.e. if you have _not_ used `input_multislice.illumination_model=4; input_multislice.temporal_spatial_incoh=1;`, which can considerably increase the simulation time and is therefore not often used), you need to "postsimulate" the effect of your probesize, by guessing on a probe shape and size, before convoluting this with your results. This can be done with

```Python
import hyperspy.api as hs
import numpy as np
from scipy.ndimage import gaussian_filter

signal = hs.load("STEM_results.hspy")

fwhm = 0.1 #FWHM [nm]
fwhm /= signal.axes_manager['x'].scale
variance = fwhm / (2 * np.sqrt( 2 * np.log( 2 ) ) ) 
blurred_signal = signal.map(gaussian_filter, inplace = False, sigma = variance, mode = 'wrap')

blurred_signal.plot()

```

### Thickness profiles
Thickness profiles can be made by creating a region of interest centered on an atomic column, and integrating the intensity inside this region. To visualise the profile of a signal with multiple detectors, you must also choose which detector(s) to show before plotting: 

```Python
import hyperspy.api as hs

signal = hs.load("STEM_results.hspy")

print(signal.metadata)
print(signal.axes_manager)

signal.plot()

roi = hs.roi.CircleROI(cx=19.07, cy=17.13, r=1.025) #Make a region of interest
print(roi)
roi.add_widget(signal, axes=['x', 'y']) #Connect the roi to the signal
cropped_signal = roi(signal) #Extract the roi from the signal (does not affect the original signal)

scattering_thickness_profile = cropped_signal.sum(axis=('x', 'y')) #Sum the signal in the x and y axes

detector = 1 #Which detector to plot the thickness profile for
scattering_thickness_profile.inav[:,detector].plot() #Plot the thickness profile for the scattering to a certain detector angle interval
```
