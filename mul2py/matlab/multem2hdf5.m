function [fid] = multem2hdf5(filename, results_struct)
%multem2hdf5 save multem data to a hdf5 file.
%   Detailed explanation goes here

%Function for writing metadata to a group.
function write_metadata(group_id, metadata, plist, filename)
        fields = fieldnames(metadata);
        for idx = 1:length(fields)
            name = fields{idx};
            value = getfield(metadata, name);
            if isequal(class(value), 'multem_input.parameters') || isequal(class(value), 'multem_input.system_conf') || isequal(class(value), 'multem_input.detector') || isequal(class(value), 'struct')
                temp_gid = H5G.create(group_id, name, plist, plist, plist);
                write_metadata(temp_gid, value, plist, filename)
                H5G.close(temp_gid);
            else
                sz = size(value);
                if sz > 0
                    acpl_id = H5P.create('H5P_ATTRIBUTE_CREATE');
                    %info = H5G.get_info(group_id);
                    if isstr(value)
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
                        %Create dataset
                        group_name = H5I.get_name(group_id);
                        if min(size(value)) < 8
                            chunks = ones(1, length(size(value)));
                        else
                            chunks = ones(1, length(size(value)))*8;
                        end                        
                        h5create(filename, sprintf('%s/%s', group_name, name), size(value), 'Datatype', 'double', 'ChunkSize', chunks, 'Deflate', 4, 'Shuffle', 1);
                        h5write(filename, sprintf('%s/%s', group_name, name), value);
%                         type_id = H5T.copy('H5T_NATIVE_DOUBLE');
%                         h5_dims = size(value);fliplr(size(value));
%                         h5_maxdims = h5_dims;%ones(size(h5_dims)) * H5ML.get_constant_value('H5S_UNLIMITED');
%                         n = size(h5_dims);
%                         n = n(2);
%                         space_id = H5S.create_simple(n, h5_dims, h5_maxdims);
%                         dcpl = 'H5P_DEFAULT';
%                         dset_id = H5D.create(group_id, name, type_id, space_id, dcpl);
%                         H5D.write(dset_id, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', plist, transpose(value));
%                         H5D.close(dset_id);
%                         H5S.close(space_id);
%                         H5T.close(type_id);
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
for idx=1:dimensions
    axis_gid = H5G.create(results_gid, sprintf('axis-%i', idx-1), plist, plist, plist);%greate group
    ax = results_struct.axes.(axes_names{dimensions - (idx-1)}); %Write the axis metadata to the group. Note that the axes_names should count "backwards"!
    write_metadata(axis_gid, ax, plist, filename);
    H5G.close(axis_gid);
end

%Create the dataset
%type_id = H5T.copy('H5T_NATIVE_DOUBLE');
%h5_dims = fliplr(size(results_struct.images));
%h5_maxdims = h5_dims;%ones(size(h5_dims)) *H5ML.get_constant_value('H5S_UNLIMITED');%* inf;
%space_id = H5S.create_simple(dimensions, h5_dims, h5_maxdims);
%dcpl = 'H5P_DATASET_CREATE';
%H5P.set_deflate(dcpl, 9);
%dset_id = H5D.create(results_gid, 'data', type_id, space_id, dcpl);
%H5D.write(dset_id, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', plist, results_struct.images);
%H5S.close(space_id);
%H5T.close(type_id);
%H5D.close(dset_id);
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
signal.metadata.signal_type = 'simulation';
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

