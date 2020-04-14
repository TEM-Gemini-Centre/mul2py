%%
clear all
clc
%% System configuration
gpu = true;
if gpu
    system_conf.precision = 1;                           % eP_Float = 1, eP_double = 2
    system_conf.device = 2;                              % eD_CPU = 1, eD_GPU = 2
    system_conf.cpu_nthread = 5; 			 % Does the number of CPU threads matter when running on GPU? EXPERIMENT!!
    system_conf.gpu_device = 1;				 % MULTEM can only use one GPU device at the time? Only ask for a single GPU from IDUN, and use this.
else
    system_conf.precision = 1;                           % eP_Float = 1, eP_double = 2
    system_conf.device = 1;                              % eD_CPU = 1, eD_GPU = 2
    system_conf.cpu_nthread = 16;
    system_conf.gpu_device = 0;
end

%% Timestamp
start_time = datetime('now','TimeZone','local');
fprintf("Starting simulation script at %s\n", start_time);

%% Paths
MULTEM_path = "/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM";  % Path to MULTEM installation on the cluster
addpath(char(sprintf("%s/crystalline_materials", MULTEM_path)));        % Add the crystalline materials function to the path (used for making models on the go if you want)
addpath(char(sprintf("%s/matlab_functions", MULTEM_path)));             % Add MULTEM matlab functions, such as `multem_default_values()` to the path.
addpath(char(sprintf("%s/mex_bin", MULTEM_path)));                      % Add the core MULTEM stuff to run simulations

%% output_details
simulation_name = "STEMcpu";
output_path = "/lustre1/work/emilc/MULTEM/Test/"; %Path to put results
mkdir(output_path);

%% Load simulation parameters. `MULTEM_input.mat` should contain a struct called `input_multislice` with all relevant simulation parameters given in its fields, including the atomistic model.
clear collection_angles
collection_angles(1).inner_ang = 0; %semi-angle (mrad)
collection_angles(1).outer_ang = 40; %semi-angle (mrad)
collection_angles(2).inner_ang = 48; %semi-angle (mrad)
collection_angles(2).outer_ang = 200; %semi-angle (mrad)

convergence_angle = 27; %semi-angle (mrad)

input_multislice = STEM_setup("Al_10x10x20.mat", convergence_angle, collection_angles, "phonons", 1, "nx", 8, "ny", 16, "instrument", "ARM200F", "multem_path", multem_path);
input_multislice.scanning_ns = 3;

%% If you want to adjust scan parameters (default is 10 pixels spanning around the center of the model with a step a/10 in x and b/10 in y:
%You can also change the scanning parameters in the call to STEM_setup()
%input_multislice.scanning_x0 = input_multislice.spec_lx/2 - input_multislice.spec_cryst_a/2;
%input_multislice.scanning_xe = input_multislice.spec_lx/2 + input_multislice.spec_cryst_a/2;

%input_multislice.scanning_y0 = input_multislice.spec_ly/2 - input_multislice.spec_cryst_b/2;
%input_multislice.scanning_ye = input_multislice.spec_ly/2 + input_multislice.spec_cryst_b/2;

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
save(sprintf("%s/%s_results.ecmat", output_path, simulation_name), "results", "-v7.3");