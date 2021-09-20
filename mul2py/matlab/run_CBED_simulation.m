function [results] = run_CBED_simulation(model_path, varargin)
    %% Timestamp
    start_time = datetime('now','TimeZone','local');
    fprintf("Starting CBED simulation function at %s\n", start_time);

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
    default_simulation_name = 'CBED';
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
    input_multem = CBED_setup(model_path, varargin{:});

    %%%%%%%%%%%%%%%%%%%%%%%%% System Setup %%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.system_conf.precision = p.Results.precision;                           % eP_Float = 1, eP_double = 2
    input_multem.system_conf.device = p.Results.device;                              % eD_CPU = 1, eD_GPU = 2
    input_multem.system_conf.cpu_nthread = p.Results.cpu_nthread;
    input_multem.system_conf.gpu_device = p.Results.gpu_device;

    %%%%%%%%%%%%%%%%%%%%%%%%% Run Simulation %%%%%%%%%%%%%%%%%%%%%%%%
    fprintf("Simulating CBED stack at (x,y) = (%f,%f)\r", input_multem.iw_x, input_multem.iw_y);
    clear il_MULTEM;
    tic;
    output_multislice = input_multem.ilc_multem;
    toc;
    end_time = datetime('now','TimeZone','local');
    fprintf("HRTEM Simulation function finished at %s\n", end_time);
    elapsed_time = seconds(end_time - start_time);
    
    %Output raw data
    parameters = input_multem.toStruct();
    save(sprintf("%s/%s_input.mat", p.Results.output_path, p.Results.simulation_name), "parameters", "-v7.3");
    save(sprintf("%s/%s_output.mat", p.Results.output_path, p.Results.simulation_name), "output_multislice", "-v7.3");
    
    %%% Create results structure %%%
    results = make_results(input_multem, output_multislice, 'elapsed_time', elapsed_time, 'title', p.Results.simulation_name);
    
    %%% Save HDF5 file %%%
    multem2hdf5(sprintf('%s/%s_results.hdf5', p.Results.output_path, p.Results.simulation_name), results);
end