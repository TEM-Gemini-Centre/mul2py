% Function to convert .mat MULTEM results to ECMAT format

function [results] = mat2ecmat(data_file, varargin)
% mat2ecmat Convert a .mat MULTEM output file to .ecmat easy-mat format for easier conversion to python with mul2py
%
%   Only works with "standard" MULTEM storage (i.e. custom structs are not likely to be treated successfully.
%
%   results = mat2ecmat('output_multislice.mat') converts the output data into a results struct with no input parameter metadata
%
%   results = mat2ecmat('output_multislice.mat', "input_parameter_file", 'input_multislice.mat') converts the output data into a results struct and adds the input parameter data from the 'input_multislice.mat' file for metadata construction in MULTEM
%
%   results = mat2ecmat('output_multislice.mat', "input_parameter_file", 'input_multislice.mat', "system_configuration_file", 'system_conf.mat') converts the data, adds the input parameters and the system configuration details for full metadata construction.

    fprintf('*** This is "mat2ecmat" converting MULTEM results to easy-mat format ***\n\n')

    p = inputParser;

    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
    validScalarNum = @(x) isnumeric(x) && isscalar(x);
    validStrChar = @(x) ischar(x) || isstring(x);
    validBool = @(x) islogical(x);
    addRequired(p, "data_file", validStrChar)
    addParameter(p, "input_parameter_file", '', validStrChar)
    addParameter(p, "system_configuration_file", '', validStrChar)
    addParameter(p, "print_details", false, validBool)
    addParameter(p, "show_progress", true, validBool)
    addParameter(p, "transpose", true, validBool)
    parse(p, data_file, varargin{:})

    [output_multislice, system_conf, input_multislice] = parse_inputs(p.Results.data_file, "input_parameter_file", p.Results.input_parameter_file, "system_configuration_file", p.Results.system_configuration_file);

    results.system = system_conf;
    results.input = input_multislice;

    try
        results.images = build_STEM_stack(output_multislice.data, "transpose", p.Results.transpose, "print_details", p.Results.print_details, "show_progress", p.Results.show_progress);
    catch ME1
        try
            results.images = build_wave_stack(output_multislice.data, "transpose", p.Results.transpose, "print_details", p.Results.print_details, "show_progress", p.Results.show_progress);
        catch ME2
            disp(ME2)
            rethrow(ME1)
        end
    end

    results.thick = output_multislice.thick;
    results.dx = output_multislice.dx;
    results.dy = output_multislice.dy;

    results.elapsed_time = nan;

    filename = strrep(p.Results.data_file, ".mat", ".ecmat");
    fprintf('\n*** Saving ecmat file to "%s" ***\n', filename)
    save(filename, "results", "-v7.3");
    fprintf('\n*** Conversion completed ***\n')

end

%Parser wrapper
function [output_multislice, system_conf, input_multislice] = parse_inputs(data_file, varargin)
    fprintf('\n*** Parsing inputs and loading files ***\n')
    p = inputParser;

    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
    validScalarNum = @(x) isnumeric(x) && isscalar(x);
    validStrChar = @(x) ischar(x) || isstring(x);

    default_system_conf = '';
    default_input_multislice = '';

    addRequired(p, "data_file", validStrChar);
    addParameter(p, "system_configuration_file", default_system_conf, validStrChar);
    addParameter(p, "input_parameter_file", default_input_multislice, validStrChar);

    parse(p, data_file, varargin{:});

    %Load output file
    fprintf('Loading data file "%s"\n', p.Results.data_file)
    data = load(p.Results.data_file);
    if length(data) == 1
        field_name = fieldnames(data);
        field_name = field_name{1};
        fprintf('\tData file "%s" only has one field with name "%s", assuming this is the data structure\n', p.Results.data_file, field_name);
        output_multislice = getfield(data, field_name);
    else
        fprintf('\tData file "%s" has multiple field names, assuming the field "output_multislice" exists\n', p.Results.data_file)
        output_multislice = data.output_multislice;
    end

    %Load system configuration
    %Default configuration
    system_conf = struct;
    system_conf.precision = nan;
    system_conf.device = nan;
    system_conf.cpu_nthread = nan;
    system_conf.gpu_device = nan;
    if ~strcmp(p.Results.system_configuration_file, default_system_conf)
        try
            fprintf('Loading system configuration "%s"\n', p.Results.system_configuration_file)
            system_conf = load(p.Results.system_configuration_file);
            if length(system_conf) == 1
                field_name = fieldnames(system_conf);
                field_name = field_name{1};
                fprintf('\tSystem configuration file "%s" only has one field with name "%s", assuming this is the configuration structure\n', p.Results.system_configuration_file, field_name);
                system_conf = getfield(system_conf, field_name);
            else
                fprintf('\tSystem configuration file "%s" has multiple field names, assuming the field "system_conf" exists\n', p.Results.system_configuration_file)
                system_conf = system_conf.system_conf;
            end
        catch ME
            fprintf('\tFailed to load system configuration\n')
            disp(ME)
        end
    end

    %Load input parameters
    %Default inputs
    input_multislice = struct;
    input_multislice.nx=nan;
    input_multislice.ny=nan;
    input_multislice.E0=nan;
    input_multislice.scanning_ns = nan;
    input_multislice.scanning_x0 = nan;
    input_multislice.scanning_y0 = nan;
    input_multislice.scanning_xe = nan;
    input_multislice.scanning_ye = nan;
    if ~strcmp(p.Results.input_parameter_file, default_input_multislice)
        try
            fprintf('Loading multislice input "%s"\n', p.Results.input_parameter_file)
            input_multislice = load(p.Results.input_parameter_file);
            if length(input_multislice) == 1
                field_name = fieldnames(input_multislice);
                field_name = field_name{1};
                fprintf('\tInput parameter file "%s" only has one field with name "%s", assuming this is the parameter structure\n', p.Results.input_parameter_file, field_name);
                input_multislice = getfield(input_multislice, field_name);
            else
                fprintf('\tInput parameter file "%s" has multiple field names, assuming the field "input_multislice" exists\n', p.Results.input_parameter_file)
                input_multislice = input_multislice.input_multislice;
            end
        catch ME
            fprintf('\tFailed to load input parameter file\n')
            disp(ME)
        end
    end
end

function print_progress(current, total, varargin)
    p = inputParser;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
    addRequired(p, "current", validScalarPosNum);
    addRequired(p, "total", validScalarPosNum);
    addParameter(p, "every", 10, validScalarPosNum);

    parse(p, current, total, varargin{:});

    if mod(p.Results.current, round(p.Results.total/p.Results.every))==0
        fprintf('\tImage %i of %i\n', p.Results.current, p.Results.total)
    end
end

%exit-wave image support functions
function [nx, ny] = get_wave_image_size(data, varargin)
    p = inputParser;
    validStruct = @(x) isstruct(x) && isfield(x, 'm2psi_tot');
    validBool = @(x) islogical(x);
    addRequired(p, "data", validStruct);
    addParameter(p, "print_details", false, validBool);
    parse(p, data, varargin{:});

    if length(data) == 1
        [nx, ny] = size(p.Results.data.m2psi_tot)
    else
        [nx, ny] = size(p.Results.data(1).m2psi_tot)
    end

    if p.Results.print_details
        fprintf('\tImage size of exit-wave data is (%ix%i) pixels\n', nx, ny);
    end
end

function image = get_wave_image(data, slice_number, varargin)
    p = inputParser;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x>0);
    validStruct = @(x) isstruct(x) && isfield(x, 'm2psi_tot');
    validBool = @(x) islogical(x);
    addRequired(p, "data", validStruct);
    addRequired(p, "slice_number", validScalarPosNum);
    addParameter(p, "transpose", true, validBool);
    addParameter(p, "print_details", false, validBool);
    parse(p, data, slice_number, varargin{:});

    try
        image = p.Results.data(p.Results.slice_number).m2psi_tot;
    catch ME
        fprintf('Could not get exit-wave image at slice number %i:\n', p.Results.slice_number)
        disp(ME)
        if slice_number==1
            fprintf('Slice number is 1. Attempting to access first slice without slice specification.\n')
            image = p.Results.data.m2psi_tot;
        else
            rethrow(ME)
        end
    end

    if p.Results.transpose
        if p.Results.print_details
            fprintf('\tTransposing image\n')
        end
        image = transpose(image);
    end
end

function image_stack = build_wave_stack(data, varargin)
    p = inputParser;
    validStruct = @(x) isstruct(x) && isfield(x, 'm2psi_tot');
    validBool = @(x) islogical(x);
    addRequired(p, "data", validStruct);
    addParameter(p, "show_progress", true, validBool);
    addParameter(p, "transpose", true, validBool);
    addParameter(p, "print_details", false, validBool);
    parse(p, data, varargin{:});

    if p.Results.print_details
        fprintf('\n*** Preparing image stack ***\n')
    end

    thicknesses = length(p.Results.data);
    if p.Results.print_details
        fprintf('\tNumber of thicknesses is %i\n', thicknesses);
    end

    %Get image size
    [nx, ny] = get_wave_image_size(p.Results.data, "print_details", p.Results.print_details);

    %Make image stack
    image_stack = zeros(nx, ny, thicknesses);
    if p.Results.print_details
        fprintf('\tImage stack size: %ix%ix%i\n', nx, ny, thicknesses)
    end

    if p.Results.print_details
        fprintf('\n *** Constructing image stack *** \n')
    end

    for t = 1:thicknesses

        if p.Results.show_progress
            print_progress(t, thicknesses)
        end

        image_stack(:, :, t) = get_wave_image(p.Results.data, t, "transpose", p.Results.transpose);
    end

    if p.Results.print_details
        fprintf('Finished building exit-wave stack\n')
    end

end

%STEM image support functions
function detectors = get_number_of_detectors(data, varargin)
    p = inputParser;
    validStruct = @(x) isstruct(x) && isfield(x, 'image_tot');
    validBool = @(x) islogical(x);
    addRequired(p, "data", validStruct);
    addParameter(p, "print_details", false, validBool);
    parse(p, data, varargin{:});

    if length(data) == 1
        detectors = length(p.Results.data.image_tot);
    else
        detectors = length(p.Results.data(1).image_tot);
    end

    if p.Results.print_details
        fprintf('Number of STEM detectors is %i\n', detectors);
    end
end

function [nx, ny] = get_STEM_image_size(data, varargin)
    p = inputParser;
    validStruct = @(x) isstruct(x) && isfield(x, 'image_tot');
    validBool = @(x) islogical(x);
    addRequired(p, "data", validStruct);
    addParameter(p, "print_details", false, validBool);
    parse(p, data, varargin{:});

    detectors = get_number_of_detectors(p.Results.data);

    if length(p.Results.data) == 1
        if detectors == 1
            [nx, ny] = size(p.Results.data.image_tot.image);
        else
            [nx, ny] = size(p.Results.data.image_tot(1).image);
        end
    else
        if detectors == 1
            [nx, ny] = size(p.Results.data(1).image_tot.image);
        else
            [nx, ny] = size(p.Results.data(1).image_tot(1).image);
        end
    end

    if p.Results.print_details
        fprintf('Image size of STEM data is (%ix%i) pixels\n', nx, ny);
    end
end

function image = get_STEM_image(data, slice_number, detector, varargin)
    p = inputParser;
    validStruct = @(x) isstruct(x) && isfield(x, 'image_tot');
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x>0);
    validBool = @(x) islogical(x);
    addRequired(p, "data", validStruct);
    addRequired(p, "slice_number", validScalarPosNum);
    addRequired(p, "detector", validScalarPosNum);
    addParameter(p, "print_details", false, validBool);
    addParameter(p, "transpose", true, validBool);
    parse(p, data, slice_number, detector, varargin{:});

    try
        image = p.Results.data(p.Results.slice_number).image_tot(p.Results.detector).image;
    catch ME
        fprintf('Could not get STEM image at slice number %i from detector %i:\n', p.Results.slice_number, p.Results.detector)
        disp(ME)
        if slice_number==1
            fprintf('Slice number is 1. Attempting to access first slice without slice specification.\n')
            image = p.Results.data.image_tot(p.Results.detector).image
        else
            rethrow(ME)
        end
    end

    if p.Results.transpose
        if p.Results.print_details
            fprintf('Transposing image\n')
        end
        image=transpose(image);
    end
end

function image_stack = build_STEM_stack(data, varargin)
    p = inputParser;
    validStruct = @(x) isstruct(x) && isfield(x, 'image_tot');
    validBool = @(x) islogical(x);
    addRequired(p, "data", validStruct);
    addParameter(p, "print_details", false, validBool);
    addParameter(p, "transpose", true, validBool);
    addParameter(p, "show_progress", true, validBool);
    parse(p, data, varargin{:});

    if p.Results.print_details
        fprintf('\n*** Preparing image stack ***\n')
    end

    thicknesses = length(p.Results.data);
    if p.Results.print_details
        fprintf('Number of thicknesses is %i\n', thicknesses);
    end

    %Get number of detectors
    detectors = get_number_of_detectors(p.Results.data, "print_details", p.Results.print_details);

    %Get image size
    [nx, ny] = get_STEM_image_size(p.Results.data, "print_details", p.Results.print_details);

    %Make image stack
    image_stack = zeros(nx, ny, thicknesses, detectors);
    if p.Results.print_details
        fprintf('Image stack size: %ix%ix%ix%i\n', nx, ny, thicknesses, detectors)
    end

    if p.Results.print_details
        fprintf('\n *** Constructing image stack *** \n')
    end

    for t = 1:thicknesses
        if p.Results.show_progress
            print_progress(t, thicknesses)
        end

        for d = 1:detectors
            image_stack(:, :, t, d) = get_STEM_image(p.Results.data, t, d, "transpose", p.Results.transpose);
        end
    end

    if p.Results.print_details
        fprintf('\n*** Finished building STEM image stack ***\n')
    end

end