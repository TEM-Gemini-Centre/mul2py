function [results] = run_STEM_simulation(model_path, alpha, collection_angles, varargin)
    %% Timestamp
    start_time = datetime('now','TimeZone','local');
    fprintf("Starting STEM simulation function at %s\n", start_time);

    %%%%%%%%%%%%%%%%%%%%%%%%% Argument Parsing %%%%%%%%%%%%%%%%%%%%%%%%
    fprintf("Parsing Arguments\n")

    %Parser object
    p = inputParser;

    %Value checkers
    validStrChar = @(x) ischar(x) || isstring(x);
    valid1or2 = @(x) x==1 || x==2;
    validPositiveNumber = @(x) isnumeric(x) && isscalar(x) && (x > 0);

    %Simulation Precision
    default_precision = 1; % eP_Float = 1, eP_double = 2
    addParameter(p, "precision", default_precision, valid1or2);

    %CPU or GPU
    default_device = 2; % eD_CPU = 1, eD_GPU = 2
    addParameter(p, "device", default_device, valid1or2)

    %Number of CPUs
    default_cpu_nthread = 1; %Number of CPUs
    addParameter(p, "cpu_nthread", default_cpu_nthread, validPositiveNumber);

    %Which GPU device to use?
    default_gpu_device = 1; % MULTEM can only use one GPU device at the time? Only ask for a single GPU from IDUN, and use this.
    addParameter(p, "gpu_device", default_gpu_device, validPositiveNumber)

    %Path to MULTEM installation
    default_MULTEM_path = '/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM';
    addParameter(p, "MULTEM_path", default_MULTEM_path, validStrChar);

    %Where to put output
    default_output_path = './'
    addParameter(p, "output_path", default_output_path, validStrChar)

    %What to label the output
    default_simulation_name = 'STEM'
    addParameter(p, "simulation_name", default_simulation_name, validStrChar)

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
    mkdir(char(output_path));

    %%%%%%%%%%%%%%%%%%%%%%%%% Simulation Setup %%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice = STEM_setup(model_path, alpha, collection_angles, varargin)

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
    n_t = length(output_multislice.data); %The number of thicknesses
    detectors = length(output_multislice.data(1).image_tot); % The number of detectors
    results.images = zeros(input_multislice.scanning_ns, input_multislice.scanning_ns, n_t, detectors);
    for t = 1:n_t
        for d = 1:detectors
            results.images(:, :, t, d) = transpose(output_multislice.data(t).image_tot(d).image);
        end
    end

    %Some additional details
    results.thick = output_multislice.thick;
    results.dx = output_multislice.dx;
    results.dy = output_multislice.dy;

    end_time = datetime('now','TimeZone','local');
    fprintf("STEM Simulation function finished at %s\n", end_time);
    results.elapsed_time = seconds(end_time - start_time);

    %% Save the data
    save(sprintf("%s/%s_results.ecmat", output_path, simulation_name), "results", "-v7.3");

end