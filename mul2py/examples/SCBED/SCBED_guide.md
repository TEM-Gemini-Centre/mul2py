# SCBED guide
This is a guide on how to prepare, setup and run a scanning convergent beam electron diffraction (SCBED) MULTEM simulation on the IDUN cluster on NTNU

## General prepartions
The first things you must consider when performing any kind of HRTEM simulation is
  - Which high-tension do you want to use?
  - Beam parameters (convergence angle)?
  - Aberrations and defocus?
  - Beam positions?
  - Simulation resolution (NOT scan-resolution)?
    - Determined by the potential sampling (number of beams)
    - Large impact on simulation time
  - When to output results?

## Making your model
After deciding on these questions, you must make your model. You should take a look at the example notebook in [https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/examples/Models](https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/examples/Models) to get a more thorough introduction to making models for MULTEM. The important thing here, is that your model size will decide (together with the potential sampling) on the maximum scattering angle you will get - and this must match your collection angles for your detectors!

## Preparing your simulation
To make the MATLAB script for running MULTEM simulations, it is probably easiest to use `mul2py`s CBED setup function: [https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/matlab/CBED_setup.m](https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/matlab/CBED_setup.m). This function will set up most of what you need, and make it easier to get going. However, if you plan on using these simulations in a publication, you should prepare your simulation yourself (based on the example scripts provided by MULTEM).

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
simulation_name = "SCBED";
output_path = "/lustre1/work/emilc/MULTEM/Test/"; %Path to put results
mkdir(output_path);
```
The first three `addpath()` statements are necessary for MATLAB to find MULTEM, and should always be included in your MULTEM scripts. The `MULTEM_path` variable in this example points to the path of the MULTEM installation on the TEM Gemini project folder on IDUN, but should be changed to your own path if you are running locally. The final three lines prepares the folder for your output and defines a name for your output files. You could also add some lines where you create and write a log file or something to test that you are allowed to write to the output folder.

### Step 3 - set up your simulation parameters
Now you can start making your `input_multislice` struct that defines the parameters used in the simulation. In this guide, we are using the `HRTEM_setup.m` function of `mul2py` to do this (you can change everything afterwards if you like):
```MATLAB
convergence_angle = 27.42 % semi-convergence angle [mrad]
input_multislice = CBED_setup("Al_10x10x20.mat", convergence_angle, "phonons", 20, "nx", 1024, "ny", 1024, "instrument", "2100F", "multem_path", multem_path);
results.input = input_multislice; % Store input parameters in results structure
results.system = system_conf;
```
You can provide several optional parameters, such as the number of phonons, the potential sampling, the instrument (predefined aberration values), and the path to the MULTEM installation, among many more.

### Step 4 - set up beam scan positions
MULTEM only supports single-position SCBED simulations. However, it is often useful to see how the beam propagates through the specimen at different positions (as a function of distance from an atomic column for instance). To do this, simply perform many MULTEM simulations in series using a for-loop. To set up this for-loop, create arrays of beam x and y positions:
```Matlab
%% Set up scan pattern
%Centre
centre_x = original_input.spec_lx/2;
centre_y = original_input.spec_ly/2;

%scan width and height
scanning_width = original_input.spec_cryst_a;
scanning_height = original_input.spec_cryst_b;

%Number of beam positions
scanning_ns_x = 2;
scanning_ns_y = 3;

%Beam start
x0 = centre_x - scanning_width/2;
y0 = centre_y - scanning_height/2;

%Beam end
xe = x0 + scanning_width;
ye = x0 + scanning_height;

%Beam positions
xs = linspace(x0, xe, scanning_ns_x);
ys = linspace(y0, ye, scanning_ns_y);

%Store beam positions in results struct
results.xs = xs;
results.ys = ys;
```

### Step 5 - Loop through beam positions
When all is prepared, you can loop through the beam positions, modify the input_multislice struct accordingly, and run the simulation by sending `system_conf` and `input_multislice` to the `il_MULTEM` function:
```matlab
results.images = zeros(input_multislice.nx, input_multislice.ny, size(xs, 2), size(ys, 2), size(input_multislice.thick, 2)); %allocate image array
results.thicknesses = {}; %cell array for storing thicknesses (in case your diffrerent beam positions have different thickensses
counter = 1; %counter to print progress
for i = 1:length(results.xs) %loop through X

    for j = 1:length(results.ys) %loop through Y

        %shift beam
        input_multislice.iw_x = xs(i);
        input_multislice.iw_y = ys(j);

        %Run simulation
        fprintf("Simulating CBED stack %i of %i: (x,y) = (%f,%f)\r", counter, length(xs) * length(ys), input_multislice.iw_x, input_multislice.iw_y); %print some progress info
        
        %Run simulation for current position
        clear il_MULTEM;
        tic;
        output_multislice = il_MULTEM(system_conf, input_multislice);
        toc;
        results.thicknesses{i, j} = output_multislice.thick; %store thicknesses
        
        %Store data
        for t = 1:length(output_multislice.data)

            results.images(:, :, i, j, t) = transpose(output_multislice.data(t).m2psi_tot); %Data must be transposed to ease import into HyperSpy

        end

        counter=counter+1; %increment counter

    end

end

%% Set some additional details
results.thick = output_multislice.thick; %Take the thicknesses of the last simulation and store them as the `thick` field for easy acess.
results.dx = output_multislice.dx;
results.dy = output_multislice.dy;

%% Save the data
save(sprintf("%s/%s_results.ecmat", output_path, simulation_name), "results", "-v7.3");
```

## Running your simulation
When you have prepared your files (the model and the matlab script), upload them to the cluster, write a suitable SLURM script, and queue the job. It is important that the SLURM script is tuned to your simulation (run time, number of cores, GPU assignment, etc) and you should make this yourself. However, the following job script should work fine for this example and serve as a template. It requires to be in the same folder as the MULTEM MATLAB script called "HRTEM.m" in this case.

```shell script
#!/bin/bash

#SBATCH --partition=GPUQ
#SBATCH --time=00-00:30:00
#SBATCH --job-name="SCBED"
#SBATCH --output=SCBED-%A.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
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

echo "Running HRTEM simulation"
matlab -nodisplay -nodesktop -nosplash -r "SCBED"

scontrol show job ${SLURM_JOB_ID} -d
```

## Converting your results
You can also use the cluster to convert your results through the `mul2py` package. The project folder on IDUN should include a virtual environment where `mul2py` is installed, and you can run the python script `convert_ecmat.py` to make a HyperSpy file from the output of the simulation.

```shell script
#!/bin/bash

#SBATCH --partition=CPUQ
#SBATCH --time=00-01:00:00
#SBATCH --job-name="SCBEDConvert"
#SBATCH --output=SCBEDConvert-%A.out
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
python /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py/mul2py/examples/convert_ecmat.py /lustre1/work/emilc/MULTEM/Test/SCBED_results.ecmat

scontrol show job ${SLURM_JOB_ID} -d
```
In this case, the shell script requires that the output from the simulation is stored at "/lustre1/work/emilc/MULTEM/Test/HRTEM_results.ecmat", and that the path to the python script "convert_ecmat.py" is correct. 

### About `convert_ecmat.py`
This script takes one argument, namely the path to a ".ecmat" file (output from the MULTEM simulation), creates a complete and scaled HyperSpy signal with sensible metadata and saves this signal to a file. The filename of the output file is the same as the input file and the extension is ".hspy". The script also takes one positional argument "--overwrite" that determines if any existing file will be overwritten (default is `True`). In other words, if you call `python convert_ecmat.py HRTEM_results.ecmat`, a file with the filename "HRTEM_results.hspy" will be created or overwritten. If you call `python convert_ecmat.py HRTEM_results.ecmat --overwrite False`, the same file will be created unless a file named "HRTEM_results.hspy" already exists in the directory.

## Postprocessing
After the output has been created (and possibly converted), it is time to make sense of it. If simply the `output_multislice` struct is written to file, it is best to do the postprocessing in MATLAB, or to export the data to some other format to treat them in your favourite tool. If an ".ecmat" file was created instead and has been converted to a ".hspy" file, you can use HyperSpy to postprocess your data.

### HyperSpy
To analyse your hyperspy SCBED simulation signal, you should first load the data, and then do whatever analysis you require. For SCBED data, you can for instance make "virtual" (in this case equivalent to actual) annular dark/bright field images:

```Python
import hyperspy.api as hs

signal = hs.load("SCBED_results.hspy") #Load the converted ".ecmat" file
print(signal.metadata) #Print the metadata of the signal
print(signal.axes_manager) #Print some info about the axes

signal.plot() #Plot the signal
roi = hs.roi.CircleROI(cx=0, cy=0, r=1.025) #Make a region of interest. Make sure that the position of the roi is within the x-y-space of the signal
roi.add_widget(signal, axes=['dx', 'dy']) #Connect the roi to the signal
cropped_signal = roi(signal) #Extract the roi from the signal (does not affect the original signal)

scattering_thickness_profile = sum(cropped_signal, axes=['dx', 'dy']) #Sum the signal in the x and y axes

scattering_thickness_profile.plot() #Plot the scattering intensity to the ROI through the thickness
```  
For a more detailed explanation and guide, please see the notebook at [https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/examples/HRTEM/SCBED_postprocessing.ipynb](https://github.com/TEM-Gemini-Centre/mul2py/tree/master/mul2py/examples/HRTEM/SCBED_postprocessing.ipynb).

## Complete example

### Model
Run the following code in e.g. a jupyter notebook to generate the model file.
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

### Simulation
Copy the following lines to a MATLAB script and save it as "HRTEM.m" in the same folder with the model you made.
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
simulation_name = "HRTEM";
output_path = "/lustre1/work/emilc/MULTEM/Test/"; %Path to put results
mkdir(output_path);

%% Load simulation parameters. `MULTEM_input.mat` should contain a struct called `input_multislice` with all relevant simulation parameters given in its fields, including the atomistic model.

input_multislice = HRTEM_setup("Al_10x10x20.mat", "phonons", 20, "nx", 1024, "ny", 1024, "instrument", "ARM200F", "multem_path", MULTEM_path);

%% Adjust defocus (set it e.g. to Scherzer defocus - this should be the default setting of `HRTEM_setup` so you must set HRTEM_setup(..., "defocus", 0) to set defocus to zero if you want no defocus.
input_multislice.obj_lens_c10 = il_scherzer_defocus(input_multislice.E_0, input_multislice.obj_lens_c_30);

%% Run simulation
clear il_MULTEM;
tic;
output_multislice = il_MULTEM(system_conf, input_multislice); 
toc;

%% Store input parameters in results struct
results.input = input_multislice;
results.system = system_conf;

%% Construct results images
results.images = zeros(input_multislice.nx, input_multislice.ny, length(output_multislice.data));
if length(output_multislice.data) == 1
    results.images(:,:, t) = transpose(output_multislice.data.m2psi_tot);
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

### SLURM
Make the following shell script and save it as "HRTEM.slurm" in the same folder as the previous two files:
```shell script
#!/bin/bash

#SBATCH --partition=GPUQ
#SBATCH --time=00-01:00:00
#SBATCH --job-name="HRTEM"
#SBATCH --output=HRTEM-%A.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
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
matlab -nodisplay -nodesktop -nosplash -r "HRTEM"

echo "Converting results to HyperSpy format using mul2py"
module load GCCcore/.8.2.0 Python/3.7.2
source /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py-env/bin/activate
python /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py/mul2py/examples/convert_ecmat.py /lustre1/work/emilc/MULTEM/Test/HRTEM_results.ecmat

scontrol show job ${SLURM_JOB_ID} -d
```
Copy both "Al_10x10x20.mat" model file and the "HRTEM.m" MATLAB script to somewhere on your directory on IDUN (see other guides on how to accomplish this). Make sure that the input path to `convert_ecmat.py` patches the path where you store the output data (specified in the MULTEM script). Then log in to IDUN on a terminal (using e.g. PuTTY on windows) and navigate to where you moved the three files "Al_10x10x20.mat", "HRTEM.m" and "HRTEM.slurm" and submit the job:
```bash
cd directory_containing_your_files
sbatch HRTEM.slurm
```
 Take a cup of coffee or maybe wait a few days for the simulation to finish before proceeding
### Conversion
Once the simulation has finished, you can convert the results on IDUN by making another SLURM script:
```shell script
#!/bin/bash

#SBATCH --partition=CPUQ
#SBATCH --time=00-01:0:00
#SBATCH --job-name="HRTEMConvert"
#SBATCH --output=HRTEMConvert-%A.out
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
python /lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py/mul2py/examples/convert_ecmat.py <directory_with_your_simulation_files>/HRTEM_results.ecmat

scontrol show job ${SLURM_JOB_ID} -d
```
Save this shell script as "HRTEMconvert.slurm" somewhere on IDUN, possibly the same place as your other simulaiton files, replace `<directory_with_your_simulation_files>` with the location of your simulation files, e.g. `/lustre1/work/emilc/MULTEM/Test/HRTEM`. Submit the job as before
 ```bash
cd <directory_with_your_simulation_files>
sbatch HRTEMconvert.slurm
```
and wait for the job to finish. Once it is finished, you can transfer the resulting `<directory_with_your_simulation_files.HRTEM_results.hspy` to your local computer and continue with postprocessing and analysing the data

### Postprocessing
Postprocess your data as you would with other TEM data, e.g. by making virtual images or thickness profile plots, e.g in a python script or a jupyter notebook.

#### Thickness profiles
Thickness profiles (i.e. extinction length profiles) can be made by creating a region of interest centered on a site, and integrating the intensity inside this region:

```Python
import hyperspy.api as hs

signal = hs.load("HRTEM_results.hspy")

print(signal.metadata)
print(signal.axes_manager)

signal.plot()

roi = hs.roi.CircleROI(cx=19.07, cy=17.13, r=1.025) #Make a region of interest
print(roi)
roi.add_widget(signal, axes=['x', 'y']) #Connect the roi to the signal
cropped_signal = roi(signal) #Extract the roi from the signal (does not affect the original signal)

scattering_thickness_profile = cropped_signal.sum(axis=('x', 'y')) #Sum the signal in the x and y axes

scattering_thickness_profile.plot() #Plot the thickness profile of a site
```
