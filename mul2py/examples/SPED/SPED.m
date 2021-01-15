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

%% Paths
MULTEM_path = "/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM";  % Path to MULTEM installation on the cluster
addpath(char(sprintf("%s/crystalline_materials", MULTEM_path)));        % Add the crystalline materials function to the path (used for making models on the go if you want)
addpath(char(sprintf("%s/matlab_functions", MULTEM_path)));             % Add MULTEM matlab functions, such as `multem_default_values()` to the path.
addpath(char(sprintf("%s/mex_bin", MULTEM_path)));                      % Add the core MULTEM stuff to run simulations

addpath(char("/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/mul2py/mul2py/matlab")) % Add mul2py matlab scripts/functions to path

%% output_details
simulation_name = "SPED";
output_path = ".";
mkdir(char(output_path));

%% Load simulation parameters. `MULTEM_input.mat` should contain a struct called `input_multislice` with all relevant simulation parameters given in its fields, including the atomistic model.
convergence_angle = 0.5; %semi-angle [mrad]
input_multislice = CBED_setup("Al_10x10x20.mat", convergence_angle, "nx", 1024, "ny", 1024, "phonons", 20, "thick_type", 2, "instrument", "2100F", "multem_path", MULTEM_path);

%% Make results struct and add inputs and system configuration
results.system = system_conf;
results.input = input_multislice;

%% Set up precession
n_theta = 360; %Cone illumination steps in degrees
phi = 1; %Cone half-angle in degrees
theta_noise_std = 0.1; #Add some randomized noise to the theta angle. Standard deviation of normalized distribution in degrees
phi_noise_std = 0.1; #Add some randomized noise to the phi angle. Standard deviation of normalized distribution in degrees

%% Set up scan pattern
centre_x = original_input.spec_lx/2;
centre_y = original_input.spec_ly/2;

scanning_width = original_input.spec_cryst_a;
scanning_height = original_input.spec_cryst_b;

scanning_ns_x = 2;
scanning_ns_y = 3;

x0 = centre_x - scanning_width/2;
y0 = centre_y - scanning_height/2;

xe = x0 + scanning_width;
ye = x0 + scanning_height;


xs = linspace(x0, xe, scanning_ns_x);
ys = linspace(y0, ye, scanning_ns_y);

%% Store scan positions in results struct
results.xs = xs;
results.ys = ys;

%% Pre-allocate image array and thickness cell
results.images = zeros(input_multislice.nx, input_multislice.ny, size(xs, 2), size(ys, 2), size(input_multislice.thick, 2));
results.thicknesses = {};

%% Loop through x and y positions
counter = 1;
for i = 1:size(results.xs, 2)
    for j = 1:size(results.ys, 2)

        %shift beam
        input_multislice.iw_x = xs(i);
        input_multislice.iw_y = ys(j);

        %Print some scan info
        fprintf("Simulating SPED stack %i of %i: (x,y) = (%f,%f)\r", counter, length(xs) * length(ys), input_multislice.iw_x, input_multislice.iw_y);
	
		%Loop through angles
		for k = 1:n_theta:
			input_multislize.theta = 360 * k / n_theta + theta_noise * randn(1);
			input_multislize.phi = phi + phi_noise * randn(1);
		
			%Run simulation
			clear il_MULTEM;
			tic;
			output_multislice = il_MULTEM(system_conf, input_multislice);
			toc;

			%Store thickness positions of output
			results.thicknesses{i, j} = output_multislice.thick;

			%Store result
			try %Try-catch to fail with some extra info and save the data so far.
				for t = 1:length(output_multislice.data)
				    results.images(:, :, i, j, t) = transpose(output_multislice.data(t).m2psi_tot); %Must be transposed to import to Hyperspy easier
				end
			catch ME
				fprintf("Exception for i=%i, j=%i, and t=%i. Data size: (%s)", i, j, t, strip(sprintf("%i,", size(output_multislice.data)), "right", ","));
				save(sprintf("%s/%s_%i_%i_%i_output.mat", output_path, simulation_name, i, j, t), "output_multislice", "-v7.3");
				save(sprintf("%s/%s_%i_%i_%i_results.ecmat", output_path, simulation_name, i, j, t), "results", "-v7.3");
				rethrow(ME)
			end
			
		end
		
    counter = counter+1;

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
