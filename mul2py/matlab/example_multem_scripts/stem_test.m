%% Configure simulation devices
system_conf.precision = 1;                           % eP_Float = 1, eP_double = 2
system_conf.device = 1;                              % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 4;                         % Number of tasks (?)
%% Load simulation parameters. `MULTEM_input.mat` should contain a struct called `input_multislice` with all relevant simulation parameters given in its fields, including the atomistic model.
clear collection_angles
collection_angles(1).inner_ang = 5; %mrad
collection_angles(1).outer_ang = 30; %mrad
collection_angles(2).inner_ang = 48; %mrad
collection_angles(2).outer_ang = 200; %mrad

convergence_angle = 15; %mrad

input_multislice = JEM2100F_STEM_setup("test_model_L_10x10x20.mat", convergence_angle, collection_angles, "phonons", 1, "nx", 256, "ny", 256, "multem_path", "C:\Program Files\MULTEM\MULTEM_binary");
input_multislice.scanning_ns = 10;

%% If you want to adjust scan parameters (default is 10 pixels spanning around the center of the model with a step a/10 in x and b/10 in y.
%input_multislice.scanning_x0 = input_multislice.spec_lx/2 - input_multislice.spec_cryst_a/2;
%input_multislice.scanning_xe = input_multislice.spec_lx/2 + input_multislice.spec_cryst_a/2;

%input_multislice.scanning_y0 = input_multislice.spec_ly/2 - input_multislice.spec_cryst_b/2;
%input_multislice.scanning_ye = input_multislice.spec_ly/2 + input_multislice.spec_cryst_b/2;

%% Store input parameters in results struct
results.input = input_multislice;

%% Run simulation
clear il_MULTEM;
tic;
output_multislice = il_MULTEM(system_conf, input_multislice); 
toc;

%% Construct results.images
n_t = length(output_multislice.data); %The number of thicknesses
detectors = length(output_multislice.data(1).image_tot); % The number of detectors

results.images = zeros(input_multislice.scanning_ns, input_multislice.scanning_ns, n_t, detectors);

try
    for t = 1:n_t
        for d = 1:detectors
            results.images(:, :, t, d) = transpose(output_multislice.data(t).image_tot(d).image);
        end
    end
catch ME
    fprintf("Exception for t=%i. Data size: (%s)", t, strip(sprintf("%i,", size(output_multislice.data)), "right", ","));
    save(sprintf("%s/%s_output.mat", results_path, result_name), "output_multislice", "-v7.3");
    save(sprintf("%s/%s_results.mat", results_path, result_name), "results", "-v7.3");
    rethrow(ME)
end

results.thick = output_multislice.thick;
results.dx = output_multislice.dx;
results.dy = output_multislice.dy;

save("stem_test_results.ecmat", "results", "-v7.3");