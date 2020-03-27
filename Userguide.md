# User guide
mul2py is made to simplify analysis of MULTEM simulation results by providing tools to convert ".mat" and ".ecmat" files to HyperSpy signals. This is very much a developing package and there are many bugs and possible pitfalls, so it should be used carefully.

In addition, mul2py provides tools for making/converting modes, making images with the atomic model overlaid the results, and to export images that can be converted to movies (in e.g. ImageJ). Several other useful tools and examples are also provided, such as scanning CBED scripts, jupyter notebooks for creating models, and sample slurm job scripts.

## Innstallation
mul2py can be found on github, and the easiest way to install it is by using git. Alternatively, the package may be downloaded manually and extracted instead of running Step 1. In addition, other environment tools than `virtualenv` may be used, such as `conda`, in which case some of the commands below must be changed. However, this is not covered at the moment.

### Step 1: Clone mul2py
Navigate to where wou want to "install" mul2py and run
```bash
git clone https://github.com/TEM-Gemini-Centre/mul2py.git
```

#### To update mul2py:
If there are updates to mul2py, this must be updated by navigating to where mul2py is installed and running
```bash
git init
git fetch origin
git checkout origin/branch
```

### Step 2: Activate/create environment
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

### Step 3: Add mul2py to environment
Add mul2py to the (active) environment by navigate to the root of mul2py and running
```bash
pip install --editable .
```
Note that the final `.` is very important.

mul2py should now be installed in editable mode and be visible when the environment is activated.

## Running MULTEM simulations
MULTEM simulations are run in MATLAB by calling
```MATLAB
output_multislice = IL_MULTEM(system_conf, input_multislice);
```
The `system_conf` and `input_multislice` are `struct` variables that are used to define the system configuration (CPU/GPU) and the input to the multislice simulation (specifying acceleration voltage, mode, specimen/atoms etc)

### Add MULTEM paths to workspace
Before you can run MULTEM, you must add the MULTEM paths to the workspace. This can be done by the following MATLAB code
```MATLAB
MULTEM_path = "/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM";  % Path to MULTEM installation on the cluster
addpath(char(sprintf("%s/crystalline_materials", MULTEM_path)));        % Add the crystalline materials function to the path (used for making models on the go if you want)
addpath(char(sprintf("%s/matlab_functions", MULTEM_path)));             % Add MULTEM matlab functions, such as `multem_default_values()` to the path.
addpath(char(sprintf("%s/mex_bin", MULTEM_path)));                      % Add the core MULTEM stuff to run simulations
```

### Specify system configuration
To run simulations, you must also tell MULTEM what kind of system configuration you want it to use (number of CPUs or if you want to use GPU). This is done by creating a `struct` and defining certain fields with specific names (these fields will be extracted later by MULTEM; so it is important that they are named correctly!). The following defines a system configuration that will use 4 CPUs.
```MATLAB
system_conf.system_conf.precision = 1;              % eP_Float = 1, eP_double = 2
system_conf.device = 1;                             % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 4;                        % number of CPUs/threads
system_conf.gpu_device = 0;                         % which GPU device to use (only if system_conf.device=2). This enables different available GPUs to be selected.
```

### Specify input parameters
Next, the input parameters to the multislice simulation must be defined. This is also done in a `struct` with specific fields. This is very sensitive to bugs, as you really have to specify the correct names. Normally, you would load the default MULTEM parameters using the `multem_default_values()` function, and then overwrite/change the fields you want/need. However, this means, that if you have typo, e.g. `input_multislice.smimulation_type=11` instead of the correct `input_multislice.simulation_type=11`, MULTEM will use whatever value that is assigned to `input_multislice.simulation_type` in `multem_default_values()`. In other words, be careful when setting up your simulation parameters!

The following code sets up a STEM simulation in MULTEM. For other simulation types, other fields in `input_Multislice` are used and must be changed instead. The best way of figuring out what parameters to use/change for the different simulation modes is by chekcing out the MULTEM examples. 

#### Set up input parameters
```MATLAB
input_multislice = multem_default_values();         % Loads default values
```

#### Set the simulation type
```MATLAB
input_multislice.simulation_type = 11; % STEM type
```

#### Set the electron-specimen interaction model
```MATLAB
input_multislice.interaction_model = 1;              % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_type = 6;                 % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6
```

#### Load a model file? See "Making models" for details
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

#### Set slicing and when to output results
```MATLAB
input_multislice.potential_slicing = 2; %Slice the potential by slice projection
input_multislice.thick_type = 2; % Get output by slices
input_multislice.thick = (0:input_multislice.spec_dz:input_multislice.spec_lz-input_multislice.spec_dz); % Get outputs for these thicknesses
```

#### Set electron-phonon interaction
```MATLAB
input_multislice.pn_model = 3;                       % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.pn_coh_contrib = 0;
input_multislice.pn_single_conf = 0;                 % 1: true, 0:false (extract single configuration)
input_multislice.pn_nconf = 20;                      % true: specific phonon configuration, false: number of frozen phonon configurations
input_multislice.pn_dim = 110;                       % phonon dimensions (xyz)
input_multislice.pn_seed = 300183;                   % Random seed(frozen phonon)
```

#### Set potential sampling.
This affects resolution/maximum scattering angle and affects the simulation time alot.
```MATLAB
input_multislice.nx = 1024;     % number of pixels in x
input_multislice.ny = 1024;     % number of pixels in y
input_multislice.bwl = 1;   % Band-width limit, 1: true, 0:false
```

#### Set acceleration voltage and beam tilt
```MATLAB
input_multislice.E_0 = 200.00;  % Acceleration Voltage (keV)
input_multislice.theta = 0.0;   % Tilt ilumination (deg)
input_multislice.phi = 0.0;     % Tilt ilumination (deg)
```

#### Set the illumination model
```MATLAB
input_multislice.illumination_model = 1;        % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
input_multislice.temporal_spatial_incoh = 1;    % 1: Temporal and Spatial, 2: Temporal, 3: Spatial
```

#### Set condenser lens parmeters
More aberrations can be defined, but C10 and C12 are the most important
```MATLAB
input_multislice.cond_lens_m = 0;                   % Vortex momentum
input_multislice.cond_lens_inner_aper_ang = 0.0;    % Inner aperture (mrad) 
input_multislice.cond_lens_outer_aper_ang = 27.42;  % Outer aperture (mrad) defines convergence semi-angle

input_multislice.obj_lens_c_10 = 0.00;  % Defocus [Å]
input_multislice.cond_lens_c12 = 0.00;  %Condenser lens spherical aberration [mm]
```

#### Set the defocus spread function (temporal incoherence?)
Only used if the illumination model is set to 4 (numerical integration) I think.
```MATLAB
dsf_sigma = il_iehwgd_2_sigma(32);                  % from defocus spread to standard deviation
input_multislice.cond_lens_dsf_sigma = dsf_sigma;   % standard deviation (�)
input_multislice.cond_lens_dsf_npoints = 5;         % # of integration points. It will be only used if illumination_model=4
```

#### Set the source spread function (spatial incoherence)
Only used if the illumination model is set to 4 (numerical integration) I think.
```MATLAB
ssf_sigma = il_hwhm_2_sigma(0.45); % half width at half maximum to standard deviation
input_multislice.cond_lens_ssf_sigma = ssf_sigma;  	% standard deviation: For parallel ilumination(�^-1); otherwise (�)
input_multislice.cond_lens_ssf_npoints = 4;         % # of integration points. It will be only used if illumination_model=4
```

#### Specify the zero defocus reference
```MATLAB
input_multislice.cond_lens_zero_defocus_type = 1;   % eZDT_First = 1, eZDT_User_Define = 2
input_multislice.cond_lens_zero_defocus_plane = 0;
```

#### Define detectors
```MATLAB
input_multislice.detector.type = 1;  % eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3

% Set collection angles (in mrad) of circular detectors. You can add as many detectors as you want by defining `input_multislice.detector.cir(x).inner_ang = <some-angle>; input_multislice.detector.cir(x).outer_ang = <some-ang>`
input_multislice.detector.cir(1).inner_ang = 5.00;      % inner collection semi-ngle of detector 1 [mrad]
input_multislice.detector.cir(1).outer_ang = 40.00;     % outer collection semi-angle of detector 1 [mrad]
input_multislice.detector.cir(2).inner_ang = 48.00;     % inner collection semi-angle of detector 2 [mrad]
input_multislice.detector.cir(2).outer_ang = 203.00;    % outer collection semi-angle of detector 2 [mrad]
```

#### Define the scan pattern
```MATLAB
input_multislice.scanning_ns = 25; % 25 probe positions
input_multislice.scanning_x0 = 0.00; % scan start in x [Å]
input_multislice.scanning_xe = 4.05; % scan stop in x [Å]
input_multislice.scanning_y0 = 0.00; % scan start in y [Å]
input_multislice.scanning_ye = 4.05; % scan stop in y [Å]
input_multislice.scanning_periodic = 0; % omit last scan row/column to enable periodic replication of images? 0: false, 1: true
```



As you can see, setting up a simulation requires quite many lines of code, and are very sensitive to typos and errors. Therefore, `mul2py` also provides some help functions to set up some usual simulation types:
 - `STEM_setup("model_path", convergence_angle, detectors)`
 - `HRTEM_setup("model_path")`
 - `CBED_setup("model_path", convergence_angle)`
 - `EWRS_setup("model_path", convergence_angle)`
 
 Each of these functions take several optional inputs to ease configuration of simulation files.
 
#### STEM_setup.m
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

#### HRTEM_setup.m
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

#### CBED_setup.m
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

#### EWRS_setup.m
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
    