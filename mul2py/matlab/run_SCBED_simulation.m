function [results] = run_SCBED_simulation(model_path, alpha, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    function [varargs] = set_vararg(varargs, name, value)
        if logical(mod(length(varargs), 2))
            error('mul2py:run_SCBED_simulation:set_vararg:OddVarargs: %s', 'Varargs should have an even number of arguments');
        end
        changed_value = false;
        for arg_idx = 1:2:length(varargs)
            if isequal(varargs{arg_idx}, name)
                varargs{arg_idx+1} = value;
                changed_value = true;
                break
            end
        end
        if ~changed_value
            varargs = [varargs {name value}];
            changed_value = true;
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

%Save separate stacs as well?
default_save_separate = 0;
addParameter(p, "save_separate", default_save_separate, validPositiveNumber);

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
    step_x = default_input.spec_cryst_a/p.Results.scan_shape(1);
else
    step_x = p.Results.step_x;
end

if p.Results.step_y == 0
    step_y = default_input.spec_cryst_b/p.Results.scan_shape(2);
else
    step_y = p.Results.step_y;
end
xs = (0:step_x:p.Results.scan_shape(1)*step_x-step_x) + default_input.iw_x; %Scan from the first scan position.
ys = (0:step_y:p.Results.scan_shape(2)*step_y-step_y) + default_input.iw_y; %Scan from the first scan position.
inputs = [];
outputs = [];

fprintf('Simulating %s using scan positions:\nx: ', p.Results.simulation_name)
fprintf('%d, ', xs);
fprintf('\ny: ')
fprintf('%d, ', ys);
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%% Run Simulation %%%%%%%%%%%%%%%%%%%%%%%%
start_time = datetime('now','TimeZone','local');
for ix = 1:p.Results.scan_shape(1)
    for iy = 1:p.Results.scan_shape(2)
        fprintf("Simulating CBED stack (%i, %i) at (x,y) = (%f,%f):\r", ix, iy, xs(ix), ys(iy));
        
        fprintf("\tUpdating input arguments...\r");
        varargin = set_vararg(varargin, 'x', xs(ix)); %update x-position of arguments
        varargin = set_vararg(varargin, 'y', ys(iy)); %update y-position of arguments
        varargin = set_vararg(varargin, 'print_details', 0); %Don't print any details.
        input = CBED_setup(model_path,  alpha, varargin{:});
        
        if p.Results.save_separate
            fprintf("\tSaving input arguments...\r");
            save(sprintf("%s/%s_input_%i_%i.mat", p.Results.output_path, p.Results.simulation_name, ix, iy), "input", "-v7.3");
        end
        
        fprintf("\tRunning simulation at (x,y) = (%f,%f)...\r", input.iw_x, input.iw_y);
        output = input.ilc_multem;
        
        if p.Results.save_separate
            fprintf("\tSaving output...\r");
            save(sprintf("%s/%s_output_%i_%i.mat", p.Results.output_path, p.Results.simulation_name, ix, iy), "output", "-v7.3");
        end
        
        fprintf("\tUpdating arrays...\r");
        inputs = [inputs input];
        outputs = [outputs output];
        
        fprintf("\tFinished with stack at (%i,%i).\r", ix, iy);
    end
end
end_time = datetime('now','TimeZone','local');
fprintf("SCBED Simulation function finished at %s\n", end_time);
elapsed_time = seconds(end_time - start_time);

fprintf("Constructing results structure\n");
%%% Create results structure %%%
results = make5Dresults(inputs, outputs, p.Results.scan_shape, xs, ys, p.Results.simulation_name, elapsed_time);

fprintf("Saving HDF5 file\n");
%%% Save HDF5 file %%%
multem2hdf5(sprintf('%s/%s_results.hdf5', p.Results.output_path, p.Results.simulation_name), results);
end
