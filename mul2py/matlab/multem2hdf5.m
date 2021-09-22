function [fid] = multem2hdf5(filename, results_struct, varargin)
%multem2hdf5 save multem results structure to a hdf5 file.
%   Writes the data and relevant metadata stored in a `struct` to a HDF5 file in a format that
%   is readable by HyperSpy (version 3.0).
%   See also:
%   `make_results`, `make5Dresults`, and `setup_axes`.

default_debug = false;
optargs = {default_debug};
numvarargs = length(varargin);
if numvarargs > length(optargs)
    error('myfuns:somefun2Alt:TooManyInputs', 'requires at most %i optional inputs', length(optargs));
end
optargs(1:numvarargs) = varargin;

[debug] = optargs{:};


    function [ext] = get_extension(filename)
        [filepath, name, ext] = fileparts(filename);
    end

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
validFilename = @(x) (ischar(x) || isstring(x)) && isequal(get_extension(x), '.hdf5');
mandatory_results_fields = ["dx" "dy" "dz" "thick" "xs" "ys" "title" "elapsed_time" "scan_shape" "input" "images" "axes"];
validResults = @(x) isstruct(x) && validate_fields(x, mandatory_results_fields);
validDebug = @(x) islogical(x);

if ~validFilename(filename)
    error('mul2py:multem2hdf5:InvalidFilename', 'Filename %s is invalid, valid filenames should contain the ".hdf5" extension.', filename);
end

if ~validResults(results_struct)
    missing_fields = [];
    results_fields = fieldnames(results_struct);
    for field_idx = 1:length(mandatory_results_fields)
        if ~ismember(mandatory_results_fields(field_idx), results_fields)
            missing_fields = [missing_fields mandandatory_results_fields(field_idx)];
        end
    end
    error('mul2py:multem2hdf5:InvalidResults', 'Provided results structure is not suited. Results structure has fields "%s" but misses fields "%s"', strjoin(fieldnames(results_struct), ', '), strjoin(missing_fields, ', '));
end

if ~validDebug(debug)
    error('mul2py:multem2hdf5:InvalidDebug', 'Debug must be a `logical`, not %s', class(debug));
end

%Function for writing metadata to a group.
    function write_metadata(group_id, metadata, plist, filename)
        if debug
            fprintf('Writing metadata (class %s) to group %s of "%s"', class(metadata), H5I.get_name(group_id), filename);
        end
        fields = fieldnames(metadata);
        for idx = 1:length(fields)
            name = fields{idx};
            value = metadata.(sprintf('%s', name));%getfield(metadata, name);
            if isequal(class(value), 'multem_input.parameters') || isequal(class(value), 'multem_input.system_conf') || isequal(class(value), 'multem_input.detector') || isequal(class(value), 'struct')
                temp_gid = H5G.create(group_id, name, plist, plist, plist);
                write_metadata(temp_gid, value, plist, filename)
                H5G.close(temp_gid);
            else
                sz = size(value);
                if sz > 0
                    group_name = H5I.get_name(group_id);
                    if debug
                        fprintf('Setting %s%s/%s of type %s\n', filename, group_name, name, class(value));
                    end
                    acpl_id = H5P.create('H5P_ATTRIBUTE_CREATE');
                    %info = H5G.get_info(group_id);
                    if ischar(value) || isstring(value)
                        type_id = H5T.copy('H5T_C_S1');
                        H5T.set_size(type_id, size(value, 2));
                        encoding = H5ML.get_constant_value('H5T_CSET_ASCII');
                        H5T.set_cset(type_id, encoding);
                        %base_id = H5T.copy('H5T_STD_I8LE');
                        %vlen_type_id = H5T.vlen_create(base_id);
                        space_id = H5S.create('H5S_SCALAR');
                        attr_id = H5A.create(group_id, name, type_id, space_id, acpl_id);
                        H5A.write(attr_id, 'H5ML_DEFAULT', value);
                        H5A.close(attr_id);
                        
                        %loc = h5info(filename, 
                        %group_name = H5I.get_name(group_id);
                        %h5writeatt(filename, group_name, name, value, 'TextEncoding', 'UTF-8')
                    elseif islogical(value)
                        parent_id = H5T.copy('H5T_NATIVE_UINT8');
                        type_id = H5T.enum_create(parent_id);
                        H5T.enum_insert(type_id, 'FALSE', 0);
                        H5T.enum_insert(type_id, 'TRUE', 1);
                        space_id = H5S.create('H5S_SCALAR');
                        attr_id = H5A.create(group_id, name, type_id, space_id, acpl_id);
                        H5A.write(attr_id, 'H5ML_DEFAULT', uint8(value));
                        H5A.close(attr_id);
                        H5S.close(space_id);
                        H5T.close(type_id);
                    elseif isscalar(value)
                        if isequal(class(value), 'int64')
                            type_id = H5T.copy('H5T_NATIVE_INT64');
                        elseif isequal(class(value), 'int32')
                            type_id = H5T.copy('H5T_NATIVE_INT32');
                        elseif isequal(class(value), 'int16')
                            type_id = H5T.copy('H5T_NATIVE_INT16');
                        elseif isequal(class(value), 'int8')
                            type_id = H5T.copy('H5T_NATIVE_INT8');
                        elseif isequal(class(value), 'uint64')
                            type_id = H5T.copy('H5T_NATIVE_UINT64');
                        elseif isequal(class(value), 'uint32')
                            type_id = H5T.copy('H5T_NATIVE_UINT32');
                        elseif isequal(class(value), 'uint16')
                            type_id = H5T.copy('H5T_NATIVE_UINT16');
                        elseif isequal(class(value), 'uint8')
                            type_id = H5T.copy('H5T_NATIVE_UINT8');
                        else
                            type_id = H5T.copy('H5T_NATIVE_DOUBLE');
                        end
                        space_id = H5S.create('H5S_SCALAR');
                        attr_id = H5A.create(group_id, name, type_id, space_id, acpl_id);
                        H5A.write(attr_id, 'H5ML_DEFAULT', value);
                        H5A.close(attr_id);
                        %group_name = H5I.get_name(group_id);
                        %h5writeatt(filename, group_name, name, value)
                    else
                        %disp(sprintf('%s: %s', name, class(value)));
                        %Create dataset
                        group_name = H5I.get_name(group_id);
                        if min(size(value)) < 8
                            chunks = ones(1, length(size(value)));
                        else
                            chunks = ones(1, length(size(value)))*8;
                        end
                        if isreal(value)
                            h5create(filename, sprintf('%s/%s', group_name, name), size(value), 'Datatype', 'double', 'ChunkSize', chunks, 'Deflate', 4, 'Shuffle', 1);
                            h5write(filename, sprintf('%s/%s', group_name, name), value);
                        else
                            h5create(filename, sprintf('%s/%s_re', group_name, name), size(real(value)), 'Datatype', 'double', 'ChunkSize', chunks, 'Deflate', 4, 'Shuffle', 1);
                            h5write(filename, sprintf('%s/%s_re', group_name, name), real(value));
                            
                            h5create(filename, sprintf('%s/%s_im', group_name, name), size(imag(value)), 'Datatype', 'double', 'ChunkSize', chunks, 'Deflate', 4, 'Shuffle', 1);
                            h5write(filename, sprintf('%s/%s_im', group_name, name), imag(value));
                        end
                    end
                    H5P.close(acpl_id);
                end
            end                
        end
    end

dimensions = size(size(results_struct.images));
dimensions = dimensions(2);


plist = 'H5P_DEFAULT';

%Create HDF5 file
fcpl = H5P.create('H5P_FILE_CREATE');
fapl = H5P.create('H5P_FILE_ACCESS');
fid = H5F.create(filename, 'H5F_ACC_TRUNC', fcpl, fapl);
%h5create(filename, dataset_location, size(results_struct.images));
%h5write(filename, dataset_location, results_struct.images);

%Set file version information
h5writeatt(filename, '/', 'file_format', 'HyperSpy');%We want to imitate the HyperSpy file format.
h5writeatt(filename, '/', 'file_format_version', '3.0'); %use 3.0 version

%Create experiments group
experiments_gid = H5G.create(fid, 'Experiments', plist, plist, plist);

%Create group for the results
results_gid = H5G.create(experiments_gid, results_struct.title, plist, plist, plist);

%Greate groups for the axes-data and set attributes
axes_names = fieldnames(results_struct.axes);
for dimension=1:dimensions
    axis_gid = H5G.create(results_gid, sprintf('axis-%i', dimension-1), plist, plist, plist);%greate group
    ax = results_struct.axes.(axes_names{dimensions - (dimension-1)}); %Write the axis metadata to the group. Note that the axes_names should count "backwards"!
    write_metadata(axis_gid, ax, plist, filename);
    H5G.close(axis_gid);
end

%Create the dataset
if results_struct.input.simulation_type == 11 || results_struct.input.simulation_type == 12
    chunks = size(results_struct.images);
else
    chunks = [32 32 ones(1, length(size(results_struct.images))-2)];
end
h5create(filename, sprintf('/Experiments/%s/data', results_struct.title), size(results_struct.images), 'Datatype', 'double', 'ChunkSize', chunks, 'Deflate', 4, 'Shuffle', 1);
h5write(filename, sprintf('/Experiments/%s/data', results_struct.title), results_struct.images);


%Create metadata groups and attributes
metadata_gid = H5G.create(results_gid, 'metadata', plist, plist, plist);
%General
general_gid = H5G.create(metadata_gid, 'General', plist, plist, plist);
general_metadata = struct('title', results_struct.title);
write_metadata(general_gid, general_metadata, plist, filename);
H5G.close(general_gid);
%Signal
signal_gid = H5G.create(metadata_gid, 'Signal', plist, plist, plist);
signal_metadata = struct();
signal_metadata.binned = 0;
signal_metadata.record_by = 'image';
signal_metadata.signal_type = 'simulation';
write_metadata(signal_gid, signal_metadata, plist, filename);
H5G.close(signal_gid);
%Parameters
%parameters_gid = H5G.create(metadata_gid, 'parameters', plist, plist, plist);
%write_metadata(parameters_gid, results_struct.input, plist, filename);
%H5G.close(parameters_gid);

%original_metadata
original_metadata_gid = H5G.create(results_gid, 'original_metadata', plist, plist, plist);
%Parameters
parameters_gid = H5G.create(original_metadata_gid, 'parameters', plist, plist, plist);
write_metadata(parameters_gid, results_struct.input, plist, filename);
if ismember('scan_inputs', fieldnames(results_struct))
    write_metadata(parameters_gid, results_struct.scan_inputs, plist, filename);
end

if ismember('probes', fieldnames(results_struct))
    if size(results_struct.probes, 1) >= 32
        x_chunks = 32;
    else
        x_chunks = size(results_struct.probes, 1);
    end
    
    if size(results_struct.probes, 2) >= 32
        y_chunks = 32;
    else
        y_chunks = size(results_struct.probes, 2);
    end
    probe_chunks = [x_chunks, y_chunks, size(results_struct.probes, 3), size(results_struct.probes, 4)];%
    loc = H5I.get_name(parameters_gid);
    h5create(filename, sprintf('%s/probes_re', loc), size(results_struct.probes), 'Datatype', 'double', 'ChunkSize', probe_chunks, 'Deflate', 4, 'Shuffle', 1);
    h5write(filename, sprintf('%s/probes_re', loc), real(results_struct.probes));
    
    h5create(filename, sprintf('%s/probes_im', loc), size(results_struct.probes), 'Datatype', 'double', 'ChunkSize', probe_chunks, 'Deflate', 4, 'Shuffle', 1);
    h5write(filename, sprintf('%s/probes_im', loc), imag(results_struct.probes));
end
H5G.close(parameters_gid);
H5G.close(original_metadata_gid);

%Close groups
H5G.close(metadata_gid);
H5G.close(results_gid);
H5G.close(experiments_gid);

%Flush data and close file.
H5F.flush(fid, 'H5F_SCOPE_GLOBAL');
H5F.close(fid);
end

