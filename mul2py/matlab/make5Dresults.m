function [results] = make5Dresults(multem_inputs, multem_outputs, scan_shape, varargin)
%MAKE5DRESULTS Make a 5D results structure based on inputs and outputs to multem.
%   make5Dresults(multem_inputs, multem_outputs, scan_shape, *xs, *ys, *title, *elapsed_time)
%   *Example*:
%
%   title = "SCBED";
%   convergence_angle = 27; %Convergence angle in mrad
%   nx = 1024; %potential sampling NOT scan points
%   ny = 1024; %potential sampling NOT scan points
%   instrument = 'ARM200F';
%   default_input = CBED_setup(model_path, convergence_angle); % for
%   calculating xs and ys
%   scan_shape = [5, 3];
%   dX = default_input.spec_cryst.a/scan_shape(1);
%   dY = default_input.spec_cryst.b/scan_shape(2);
%   xs = [0:dX:scan_shape(1)*dX] + default_input.iw_x; %Scan from the first scan position
%   ys = [0:dY:scan_shape(2)*dY] + default_input.iw_y; %Scan from the first scan position.
%   inputs = [];
%   outputs = [];
%   start_time = datetime('now','TimeZone','local');
%   for ix = 1:scan_shape(1)
%       for iy = 1:scan_shape(2)
%           input = CBED_setup(model_path, convergence_angle, 'x', xs(ix), 'y', ys(iy), 'nx', nx, 'ny', ny, 'instrument', instrument);
%           save(sprintf("%s_input_%i_%i.mat", title, ix, iy), "input", "-v7.3");
%           output = input.ilc_multem;
%           save(sprintf("%s_output_%i_%i.mat", title, ix, iy), "output", "-v7.3");
%           inputs = [inputs input];
%           outputs = [outputs output];
%   end_time = datetime('now','TimeZone','local');
%   results = make5Dresults(inputs, outputs, scan_shape, xs, ys, "SCBED", seconds(end_time - start_time))


% Validation functions
    function [valid] = validate_fields(structure, mandatory_fields)
        structure_fields = fieldnames(structure);
        for name_idx = 1:length(mandatory_fields)
            valid = ismember(mandatory_fields(name_idx), structure_fields);
            if ~valid
                return
            end
        end
        valid = true;
    end

    function [valid] = validate_outputs(output_structure)
        mandatory_output_fields = ["dx", "dy", "x", "y", "thick", "data"];
        valid = validate_fields(output_structure, mandatory_output_fields);
        if ~valid
            return
        end
        for idx=1:length(output_structure)
            valid = ismember('m2psi_tot', fieldnames(output_structure(idx).data));
            if ~valid
                return
            end
        end
    end
    
    function [struct_array] = reshape_struct_vector(struct_vector, shape)
        if ~isequal(prod(shape), length(struct_vector))
            error('mul2py:run_SCBED_simulation:reshape_struct_vector:DimnsionError: %s', 'prod(%f) of shape does not match length %i of struct_vector', prod(shape), length(struct_vector));
        end

        struct_array = struct_vector(1);

        counter = 1;
        for x = 1:shape(1)
            for y = 1:shape(2)
                struct_array(x, y) = struct_vector(counter);
                counter = counter + 1 ;
            end
        end
    end

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
validStrChar = @(x) ischar(x) || isstring(x);
validVector = @(x) isrow(x) || iscol(x);
mandatory_input_fields = fieldnames(multem_input.parameters);
validMultemInputs = @(x) length(x) >= 1 && (isequal(class(x), 'multem_input.parameters') || isstruct(x)) && validate_fields(x, mandatory_input_fields);
validMultemOutputs = @(x) length(x) >= 1 && isstruct(x) && validate_outputs(x);
validShape = @(x) isnumeric(x) && (isrow(x) || iscolumn(x)) && length(x)==2;


%%%%%%%%%%%%%%%%%%%%%%%%% Argument Parsing %%%%%%%%%%%%%%%%%%%%%%%%

%Validate required inputs
if ~validMultemInputs(multem_inputs)
    error('mul2py:multem2hdf5:InvalidResults', 'Provided results structure of class "%s" with length %i is not supported', class(multem_inputs), length(multem_inputs));
end
if ~validMultemOutputs(multem_outputs)
    error('mul2py:multem2hdf5:InvalidOutputs', 'Provided results outputs of class "%s" is not supported', class(multem_outputs));
end
if ~validShape(scan_shape)
    error('mul2py:multem2hdf5:InvalidScanShape', 'Provided scan shape of class "%s" and length %i is not supported', class(scan_shape), length(scan_shape));
end

%validate optional inputs
default_xs = 1:scan_shape(1);
default_ys = 1:scan_shape(2);
default_title = 'MULTEM_5Dsimulation';
default_elapsed_time = 0;

%Construct optional arguments
optargs = {default_xs default_ys default_title default_elapsed_time};
numvarargs = length(varargin);
if numvarargs > length(optargs)
    error('mul2py:make5Dresults:TooManyInputs', 'requires at most %i optional inputs', length(optargs));
end
optargs(1:numvarargs) = varargin;

[xs, ys, title, elapsed_time] = optargs{:};

if ~(validVector(xs) && length(xs) == scan_shape(1))
    error('mul2py:multem2hdf5:InvalidScanpositions', 'Provided x scan positions of class "%s" and length %i is not supported for scan shape along x equals to %i', class(xs), length(xs), scan_shape(1));
end
if ~validVector(ys)
    error('mul2py:multem2hdf5:InvalidScanpositions', 'Provided y scan positions of class "%s" and length %i is not supported for scan shape along x equals to %i', class(ys), length(ys), scan_shape(2));
end
if ~(length(xs)*length(ys) == scan_shape(1)*scan_shape(2))
    error('mul2py:multem2hdf5:IncompatibleScanDimensions', 'Provided scan positions of length %ix%i does not match provided scan shape %ix%i', length(xs), length(ys), scan_shape(1), scan_shape(2));
end
if ~validStrChar(title)
    error('mul2py:multem2hdf5:InvalidTitle', 'Provided title of class "%s" is invalid. Expected a string or char.', class(title));
end
if ~validScalarPosNum(elapsed_time)
    error('mul2py:multem2hdf5:InvalidElapsedTime', 'Provided elapsed time of class "%s" is invalid. Expected a positive number', class(elapsed_time));
end

%Create results
results = struct();

%Store metadata
%results.simulation_type = p.Results.multem_inputs(1).simulation_type;
results.dz = multem_inputs(1).spec_dz;
results.thick = multem_outputs(1).thick;

%Find scan positions
if isnan(xs)
    results.xs = 1:scan_shape(1);
elseif length(xs) ~= scan_shape(1)
    if length(xs) == 1
        dX = xs;
    else
        dX = xs(2) - xs(1);
    end
    xmax = dX*scan_shape(1);
    results.xs = dX:dX:xmax;
else
    results.xs = xs;
end

if isnan(ys)
    results.ys = 1:scan_shape(2);
elseif length(ys) ~= scan_shape(2)
    if length(ys) == 1
        dY = ys;
    else
        dY = ys(2) - ys(1);
    end
    ymax = dY*scan_shape(2);
    results.ys = dY:dY:ymax;
else
    results.ys = ys;
end

results.title = title;
results.elapsed_time = elapsed_time;
results.dx = multem_outputs(1).dx;
results.dy = multem_outputs(1).dy;
results.scan_shape = scan_shape;

%Store inputs of first simulation only. THe other inputs *should* be very
%similar (except for iw_psi_x/y).
results.input = multem_inputs(1);
for i=1:length(multem_inputs)
    input = multem_inputs(i);
    if isequal(class(multem_inputs(i)), 'multem_input.parameters')
        input = input.toStruct();
    end
    results.scan_inputs.(sprintf('simulation%i', i)) = input;
end

%Create array for storing probes
results.probes = zeros(multem_inputs(1).nx, multem_inputs(1).ny, scan_shape(1), scan_shape(2));

%Construct data array
results.images = zeros(multem_inputs(1).nx, multem_inputs(1).ny, scan_shape(1), scan_shape(2), length(results.thick));
outputs = reshape_struct_vector(multem_outputs, scan_shape);
inputs = reshape_struct_vector(multem_inputs, scan_shape);
for x = 1:scan_shape(1)
    for y = 1:scan_shape(2)
        output = outputs(x, y);
        if length(output.data) == 1
            results.images(:,:, 1, 1, 1) = transpose(output.data.m2psi_tot);
        else
            for t = 1:length(output.data)
                results.images(:, :, x, y, t) = transpose(output.data(t).m2psi_tot);
            end
        end
        results.probes(:, :, x, y) = transpose(inputs(x, y).iw_psi);
    end
end

results.axes = setup_axes(results);
end

