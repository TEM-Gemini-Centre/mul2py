function [results] = run_HRTEM_simulation(model_path, varargin)
    %% Timestamp
    start_time = datetime('now','TimeZone','local');
    fprintf("Starting HRTEM simulation function at %s\n", start_time);

    %%%%%%%%%%%%%%%%%%%%%%%%% Argument Parsing %%%%%%%%%%%%%%%%%%%%%%%%
    fprintf("Parsing Arguments\n")

    %Parser object
    p = inputParser;
    p.KeepUnmatched = true;

    %Value checkers
    validStrChar = @(x) ischar(x) || isstring(x);
    valid1or2 = @(x) x==1 || x==2;
    validPositiveNumber = @(x) isnumeric(x) && isscalar(x) && (x > 0);

    %Simulation Precision
    default_precision = 1; % eP_Float = 1, eP_double = 2
    addParameter(p, "precision", default_precision, valid1or2);

    %CPU or GPU
    default_device = 2; % eD_CPU = 1, eD_GPU = 2
    addParameter(p, "device", default_device, valid1or2);

    %Number of CPUs
    default_cpu_nthread = 1; %Number of CPUs
    addParameter(p, "cpu_nthread", default_cpu_nthread, validPositiveNumber);

    %Which GPU device to use?
    default_gpu_device = 1; % MULTEM can only use one GPU device at the time? Only ask for a single GPU from IDUN, and use this.
    addParameter(p, "gpu_device", default_gpu_device, validPositiveNumber);

    %Path to MULTEM installation
    default_MULTEM_path = '/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM';
    addParameter(p, "MULTEM_path", default_MULTEM_path, validStrChar);

    %Save results, or only return them?
    default_save = 1;
    addParameter(p, "save", default_save, validPositiveNumber);

    %Where to put output
    default_output_path = './';
    addParameter(p, "output_path", default_output_path, validStrChar);

    %What to label the output
    default_simulation_name = 'HRTEM';
    addParameter(p, "simulation_name", default_simulation_name, validStrChar);

    %Print parser info?
    default_print_parser = 0;
    addParameter(p, "print_parser", default_print_parser, validPositiveNumber);

    %Parse the arguments
    parse(p, varargin{:});

    if p.Results.print_parser
        fprintf("Parser:\n")
        disp(p)

        fprintf("Parser results:\n")
        disp(p.Results)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%% Add MULTEM paths %%%%%%%%%%%%%%%%%%%%%%%%
    addpath(char(sprintf("%s/crystalline_materials", p.Results.MULTEM_path)));
    addpath(char(sprintf("%s/matlab_functions", p.Results.MULTEM_path)));
    addpath(char(sprintf("%s/mex_bin", p.Results.MULTEM_path)));

    %%%%%%%%%%%%%%%%%%%%%%%%% Make Output Path %%%%%%%%%%%%%%%%%%%%%%%%
    if p.Results.save
        mkdir(char(p.Results.output_path));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%% Simulation Setup %%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice = HRTEM_setup(model_path, varargin{:});

    %%%%%%%%%%%%%%%%%%%%%%%%% System Setup %%%%%%%%%%%%%%%%%%%%%%%%
    system_conf.precision = p.Results.precision;                           % eP_Float = 1, eP_double = 2
    system_conf.device = p.Results.device;                              % eD_CPU = 1, eD_GPU = 2
    system_conf.cpu_nthread = p.Results.cpu_nthread;
    system_conf.gpu_device = p.Results.gpu_device;

    %%%%%%%%%%%%%%%%%%%%%%%%% Run Simulation %%%%%%%%%%%%%%%%%%%%%%%%
    fprintf("Running simulation\n")
    clear il_MULTEM;
    tic;
    output_multislice = il_MULTEM(system_conf, input_multislice);
    toc;

    %%%%%%%%%%%%%%%%%%%%%%%%% Construct results %%%%%%%%%%%%%%%%%%%%%%%%
    fprintf("Constructing results structure\n")
    %Input Parameters
    results.input = input_multislice;
    results.system = system_conf;

    %Images
    results.images = zeros(input_multislice.nx, input_multislice.ny, length(output_multislice.data));
    if length(output_multislice.data) == 1
        results.images(:,:, 1) = transpose(output_multislice.data.m2psi_tot);
    else
        for t = 1:length(output_multislice.data)
            results.images(:, :, t) = transpose(output_multislice.data(t).m2psi_tot);
        end
    end

    %Some additional details
    results.thick = output_multislice.thick;
    results.dx = output_multislice.dx;
    results.dy = output_multislice.dy;

    end_time = datetime('now','TimeZone','local');
    fprintf("HRTEM Simulation function finished at %s\n", end_time);
    results.elapsed_time = seconds(end_time - start_time);

    %% Save the data
    if p.Results.save
        save(sprintf("%s/%s_results.ecmat", p.Results.output_path, p.Results.simulation_name), "results", "-v7.3");
    end

end