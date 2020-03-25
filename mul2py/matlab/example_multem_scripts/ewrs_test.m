%% Configure simulation devices
system_conf.precision = 1;                           % eP_Float = 1, eP_double = 2
system_conf.device = 1;                              % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 4;                         % Number of tasks (?)
%% Load simulation parameters. `MULTEM_input.mat` should contain a struct called `input_multislice` with all relevant simulation parameters given in its fields, including the atomistic model.
input_multislice = JEM2100F_EWRS_setup("test_model_L_10x10x20.mat", 27.42, "nx", 256, "ny", 256, "multem_path", "C:\Program Files\MULTEM\MULTEM_binary");
original_input = input_multislice;
%% Set up scan pattern
%Scan one unit cell centered in the middle of the cell with 1/4 unit cell resolution
centre_x = original_input.spec_lx/2;
centre_y = original_input.spec_ly/2;

step_x = original_input.spec_cryst_a/4;
step_y = original_input.spec_cryst_b/4;

width = 3; %Steps
height = 3; %Steps

xs = (centre_x - width / 2 * step_x : step_x : centre_x + width / 2 * step_x);
ys = (centre_y - height / 2 * step_y : step_y : centre_y + height / 2 * step_y);

%% Loop through x and y positions
results.input = original_input;
results.xs = xs;
results.ys = ys;
results.images = zeros(input_multislice.nx, input_multislice.ny, size(xs, 2), size(ys, 2), size(input_multislice.thick, 2));
results.thicknesses = {};
counter = 1;
for i = 1:size(results.xs, 2)
    x = xs(i);
    for j = 1:size(results.ys, 2)
        y = ys(j);
        %shift beam
        input_multislice = original_input;
        input_multislice.iw_x = x;
        input_multislice.iw_y = y;
        
        file_name = sprintf("%s_%i_%i", result_name, i, j); %Name for storing individual output file if an exception is thrown.

        %Run simulation
        clear il_MULTEM;
        fprintf("Simulating exit wave (real-space) stack %i of %i: (x,y) = (%f,%f)\r", counter, size(xs, 2)*size(ys, 2), input_multislice.iw_x, input_multislice.iw_y);
        tic;
        output_multislice = il_MULTEM(system_conf, input_multislice); 
        toc;
        results.thicknesses{i, j} = output_multislice.thick;
        try
            for t = 1:size(output_multislice.data, 2)
                results.images(:, :, i, j, t) = transpose(output_multislice.data(t).m2psi_tot);
            end
        catch ME
            fprintf("Exception for i=%i, j=%i, and t=%i. Data size: (%s)", i, j, t, strip(sprintf("%i,", size(output_multislice.data)), "right", ","));
            save(sprintf("%s/%s_output.mat", results_path, file_name), "output_multislice", "-v7.3");
            save(sprintf("%s/%s_results.mat", results_path, result_name), "results", "-v7.3");
            rethrow(ME)
        end        
        counter=counter+1;
    end
end
results.dx = output_multislice.dx;
results.dy = output_multislice.dy;

save("ewrs_test_results.ecmat", "results", "-v7.3");