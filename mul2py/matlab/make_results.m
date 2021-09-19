function [results] = make_results(multem_input, multem_output, varargin)
%make_results Make a results structure based on input and output to multem.
%   Constructs a struct with detailed metadata suited for writing to HDF5
%   using `multem2hdf5.m`. When used together with custom scanning schemes,
%   `xs` and `ys` should be specified as this will generate the proper size
%   of the image-array. In such cases, this function should only be called
%   for the first output, and subsequent outputs should be stored manually
%   in the proper location of the data array.


%%%%%%%%%%%%%%%%%%%%%%%%% Argument Parsing %%%%%%%%%%%%%%%%%%%%%%%%
default_xs = 1;
default_ys = 1;
default_title = 'MULTEM_simulation';
default_elapsed_time = 0;

p = inputParser;
p.KeepUnmatched = true;

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
validStrChar = @(x) ischar(x) || isstring(x);
validVector = @(x) isrow(x) || iscol(x);
validMultemInput = @(x) isequal(class(x), 'multem_input.parameters') || isstruct(x);
validMultemOutput = @(x) isstruct(x) && ismember('data', fieldnames(x)) && ismember('dx', fieldnames(x)) && ismember('dy', fieldnames(x)) && (ismember('image_tot', fieldnames(x.data)) || ismember('m2psi_tot', fieldnames(x.data)));


addRequired(p, "multem_input", validMultemInput);
addRequired(p, "multem_output", validMultemOutput);
addParameter(p, "xs", default_xs, validVector);
addParameter(p, "ys", default_ys, validVector);
addParameter(p, "title", default_title, validStrChar);
addParameter(p, "elapsed_time", default_elapsed_time, validScalarPosNum);
parse(p, multem_input, multem_output, varargin{:});
    
    
results = struct();
if isequal(class(multem_input, 'multem_input.parameters'))
    results.input = multem_input.toStruct();
else
    results.input = multem_input;
end
results.dz = multem_input.spec_dz;
results.thick = multem_output.thick;
results.xs = p.Results.xs;
results.ys = p.Results.ys;
results.title = p.Results.title;
results.elapsed_time = p.Results.elapsed_time;

if isequal(multem_input.simulation_type, 11) || isequal(multem_input.simulation_type, 12)
    %STEM
    %Find number of detectors
    if input_multem.detector.type == 1
        detectors = length(multem_input.detector.cir);
    elseif input_multem.detector.type == 2
        detectors = length(multem_input.detector.radial);
    elseif input_multem.detector.type == 3
        detectors = length(multem_input.detector.matrix);
    else
        detectors = length(multem_output.data(1).image_tot);
    end
    
    results.dx = (input_multem.scanning.xe-input_multem.scanning.x0) / input_multem.scanning_ns;
    results.dy = (input_multem.scanning.ye-input_multem.scanning.y0) / input_multem.scanning_ns;
    results.images = zeros(input_multem.scanning_ns, input_multem.scanning_ns, length(multem_output.data), detectors);
    if length(multem_output.data) == 1
        for det=1:detectors
            results.images(:,:, 1, det) = transpose(multem_output.data.image_tot(det).image);
        end
    else
        for t = 1:length(multem_output.data)
            if length(multem_output.data(t).image_tot) == 1
                results.images(:, :, t, 1) = transpose(multem_output.data(t).image_tot.image);
            else
                for det=1:detectors
                    results.images(:, :, t, det) = transpose(multem_output.data(t).image_tot(det).image);
                end
            end
        end
    end
else
    %HRTEM and other direct data: Should not be used for custom "scanning"
    %modes
    results.dx = multem_output.dx;
    results.dy = multem_output.dy;
    if length(p.Results.xs) > 1 && length(p.Results.ys) > 1
        results.images = zeros(multem_input.nx, multem_input.ny, length(multem_output.data), length(p.Results.xs), length(p.Results.ys));
        if length(multem_output.data) == 1
            results.images(:,:, 1, 1, 1) = transpose(multem_output.data.m2psi_tot);
        else
            for t = 1:length(multem_output.data)
                results.images(:, :, t, 1, 1) = transpose(multem_output.data(t).m2psi_tot);
            end
        end
    elseif length(p.Results.xs) > 1
        results.images = zeros(multem_input.nx, multem_input.ny, length(multem_output.data), length(p.Results.xs));
        if length(multem_output.data) == 1
            results.images(:,:, 1, 1) = transpose(multem_output.data.m2psi_tot);
        else
            for t = 1:length(multem_output.data)
                results.images(:, :, t, 1) = transpose(multem_output.data(t).m2psi_tot);
            end
        end
    elseif length(p.Results.ys) > 1
        results.images = zeros(multem_input.nx, multem_input.ny, length(multem_output.data), length(p.Results.ys));
        if length(multem_output.data) == 1
            results.images(:,:, 1, 1) = transpose(multem_output.data.m2psi_tot);
        else
            for t = 1:length(multem_output.data)
                results.images(:, :, t, 1) = transpose(multem_output.data(t).m2psi_tot);
            end
        end
    else
        results.images = zeros(multem_input.nx, multem_input.ny, length(multem_output.data));
        if length(multem_output.data) == 1
            results.images(:,:, 1) = transpose(multem_output.data.m2psi_tot);
        else
            for t = 1:length(multem_output.data)
                results.images(:, :, t) = transpose(multem_output.data(t).m2psi_tot);
            end
        end
    end
end

results.axes = setup_axes(results);
    
end

