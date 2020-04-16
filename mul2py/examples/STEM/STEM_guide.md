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
```MULTEM
input_multislice.scanning_ns = 25;
input_multislice.scanning_periodic = 0;
input_multislice.scanning_x0 = input_multislice.spec_lx/2 - input_multislice.spec_cryst_a/2;
input_multislice.scanning_xe = input_multislice.spec_lx/2 + input_multislice.spec_cryst_a/2;

input_multislice.scanning_y0 = input_multislice.spec_ly/2 - input_multislice.spec_cryst_b/2;
input_multislice.scanning_ye = input_multislice.spec_ly/2 + input_multislice.spec_cryst_b/2;
```
This will scan 25 pixels (`input_multislice.scanning_ns=25`) and include the last row and column (`input_multislice.scanning_periodic=0`). It will start half a unit cell from the center of the model (`input_multislice.scanning_x0 = input_multislice.spec_lx/2 - input_multislice.spec_cryst_a/2`) and scan one unit cell to end up at half a unit cell away from the centre of the model in the other direction (`input_multislice.scanning_x0 = input_multislice.spec_lx/2 + input_multislice.spec_cryst_a/2`).

### Step 5 - Run simulation
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

## Postprocessing

## Complete example
