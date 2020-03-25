%% Configure simulation devices
system_conf.precision = 1;                           % eP_Float = 1, eP_double = 2
system_conf.device = 1;                              % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 4;                         % Number of tasks (?)
%% Load simulation parameters. `MULTEM_input.mat` should contain a struct called `input_multislice` with all relevant simulation parameters given in its fields, including the atomistic model.
input_multislice = JEM2100F_CBED_setup("test_model_L_10x10x20.mat", 7.5, "nx", 256, "ny", 256, "multem_path", "C:\Program Files\MULTEM\MULTEM_binary");
original_input = input_multislice;

%% Create results struct
results.input = input_multislice;
results.x = input_multislice.iw_x;
results.y = input_multislice.iw_y;
%%
clear il_MULTEM;
tic;
output_multislice = il_MULTEM(system_conf, input_multislice); 
toc;

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
results.thicknesses = output_multislice.thick;

results.dx = output_multislice.dx;
results.dy = output_multislice.dy;
%%
save"cbed_test_results.ecmat", "results", "-v7.3");
%%