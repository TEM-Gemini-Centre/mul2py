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
    system_conf.cpu_nthread = 4; 
    system_conf.gpu_device = 0;
end

%% Timestamp
start_time = datetime('now','TimeZone','local');
fprintf("Starting simulation script at %s\n", start_time);

%% output_details
simulation_name = "CBED_test";
output_path = ".";
mkdir(output_path);

%% Load simulation parameters. `MULTEM_input.mat` should contain a struct called `input_multislice` with all relevant simulation parameters given in its fields, including the atomistic model.
input_multislice = CBED_setup("test_model_L_10x10x20.mat", 7.5, "phonons", 1, "nx", 8, "ny", 16, "instrument", "2100F", "multem_path", "C:\Program Files\MULTEM\MULTEM_binary");
original_input = input_multislice;

%% Create results struct
results.input = input_multislice;
results.x = input_multislice.iw_x;
results.y = input_multislice.iw_y;

%% Run simulation
clear il_MULTEM;
tic;
output_multislice = il_MULTEM(system_conf, input_multislice); 
toc;

%% Construct results.images

results.images = zeros(input_multislice.nx, input_multislice.ny, length(output_multislice.thick));

if length(output_multislice.data) == 1
    results.images(:,:, 1) = transpose(output_multislice.data.m2psi_tot);
else
    try
        for t = 1:length(output_multislice.data)
            results.images(:, :, t) = transpose(output_multislice.data(t).m2psi_tot);
        end
    catch ME
        fprintf("Exception for t=%i. Data size: (%s)", t, strip(sprintf("%i,", size(output_multislice.data)), "right", ","));
        save(sprintf("%s/%s_output.mat", results_path, result_name), "output_multislice", "-v7.3");
        save(sprintf("%s/%s_results.mat", results_path, result_name), "results", "-v7.3");
        rethrow(ME)
    end
end

%% Set some additional details
results.thick = output_multislice.thick;
results.dx = output_multislice.dx;
results.dy = output_multislice.dy;

results.simulation_type = "CBED"; %This will be used by mul2py to decide how to construct HyperSpy signals

end_time = datetime('now','TimeZone','local');
fprintf("Simulation finished at %s\n", end_time);
results.elapsed_time = seconds(end_time - start_time);

%% Save the data
save(sprintf("%s/%s_results.ecmat", output_path, simulation_name), "results", "-v7.3");