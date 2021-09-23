function [results] = run_SCBED_simulation(model_path, alpha, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    function [varargs] = set_vararg(varargs, name, value)
        if logical(mod(length(varargs), 2))
            error('mul2py:run_SCBED_simulation:set_vararg:OddVarargs: %s', 'Varargs should have an even number of arguments');
        end
        arg_idx = 1;
        while arg_idx <= length(varargs)
            if isequal(varargs{arg_idx}, name)
                varargs{arg_idx+1} = value;
                break
            end
        end
    end
%% Timestamp
start_time = datetime('now','TimeZone','local');
fprintf("Starting SCBED simulation function at %s\n", start_time);

%%%%%%%%%%%%%%%%%%%%%%%%% Argument Parsing %%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Parsing Arguments\n")

%Parser object
p = inputParser;
p.KeepUnmatched = true;

%Value checkers
validStrChar = @(x) ischar(x) || isstring(x);
validPositiveNumber = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
validPositiveNumberOrNan = @(x) isnumeric(x) && isscalar(x) && ((x >= 0) || isnan(x));
validShape = @(x) isnumeric(x) && (isrow(x) || iscolumn(x)) && length(x)==2;

%x-scan
default_step_x = 0;
addParameter(p, 'step_x', default_step_x, validPositiveNumber);

%y-scan
default_step_y = 0;
addParameter(p, 'step_y', default_step_y, validPositiveNumber);

%scan_shape
default_scan_shape = [3 3];
addParameter(p, 'scan_shape', default_scan_shape, validShape);

%scan x start
default_x0 = nan;
addParameter(p, 'x0', default_x0, validPositiveNumberOrNan);

%scan y start
default_y0 = nan;
addParameter(p, 'y0', default_y0, validPositiveNumberOrNan);

%Path to MULTEM installation
default_MULTEM_path = '/cluster/projects/itea_lille-nv-fys-tem/repositories/multem';
addParameter(p, "MULTEM_path", default_MULTEM_path, validStrChar);

%Save results, or only return them?
default_save = 1;
addParameter(p, "save", default_save, validPositiveNumber);

%Where to put output
default_output_path = './';
addParameter(p, "output_path", default_output_path, validStrChar);

%What to label the output
default_simulation_name = 'SCBED';
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
default_input = CBED_setup(model_path, alpha, varargin{:}); %Get default input to calculate scan region if default values
if p.Results.step_x == 0
    step_x = default_input.spec_cryst.a/p.Results.scan_shape(1);
else
    step_x = p.Results.step_x;
end

if p.Results.step_y == 0
    step_y = default_input.spec_cryst.b/p.Results.scan_shape(2);
else
    step_y = p.Results.step_y;
end
xs = (0:step_x:p.Results.scan_shape(1)*step_x) + default_input.iw_x; %Scan from the first scan position.
ys = (0:step_y:p.Results.scan_shape(2)*step_y) + default_input.iw_y; %Scan from the first scan position.
inputs = [];
outputs = [];

fprintf('Simulating %s using scan positions:\nx: ', p.Results.title)
fprintf('%d, ', xs);
fprintf('\ny: ')
fprintf('%d, ', ys);
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%% Run Simulation %%%%%%%%%%%%%%%%%%%%%%%%
start_time = datetime('now','TimeZone','local');
for ix = 1:p.Results.scan_shape(1)
    for iy = 1:p.Results.scan_shape(2)
        varargin = set_vararg(varargin{:}, 'x', xs(ix)); %update x-position of arguments
        varargin = set_vararg(varargin{:}, 'y', ys(iy)); %update y-position of arguments
        input = CBED_setup(model_path,  alpha, varargin{:});
        save(sprintf("%s_input_%i_%i.mat", title, ix, iy), "input", "-v7.3");
        fprintf("Simulating CBED stack at (x,y) = (%f,%f)\r", input_multem.iw_x, input_multem.iw_y);
        output = input.ilc_multem;
        save(sprintf("%s_output_%i_%i.mat", title, ix, iy), "output", "-v7.3");
        inputs = [inputs input];
        outputs = [outputs output];
    end
end
end_time = datetime('now','TimeZone','local');
fprintf("SCBED Simulation function finished at %s\n", end_time);
elapsed_time = seconds(end_time - start_time);

%%% Create results structure %%%
results = make5Dresults(inputs, outputs, scan_shape, xs, ys, p.Results.simulation_name, elapsed_time);

%%% Save HDF5 file %%%
multem2hdf5(sprintf('%s/%s_results.hdf5', p.Results.output_path, p.Results.simulation_name), results);
end

