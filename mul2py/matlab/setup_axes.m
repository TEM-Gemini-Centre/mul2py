function [axes_struct] = setup_axes(results_struct)
%setup_axes Setup the axes for the provided multem simulation results
%structure.
%   The results structure should have the data array stored in the "images"
%   field, the data scales (scale of the simulation data) should be stored
%   in `dx` and `dy` fields, the thicknesses in the `thick` field, and the
%   input parameters in the `input` field.

axes_struct = struct();

dimensions = size(size(results_struct.images), 2);

%Determine predefined dimension catagories
if isequal(results_struct.input.simulation_type, 11) || isequal(results_struct.input.simulation_type, 12)
    if dimensions >=5 %5D dataset
        data_x_dim = 1;
        data_y_dim = 2;
        data_X_dim = 3;
        data_Y_dim = 4;
        data_z_dim = 5;
        data_detector_dim = 6;
    else
        data_x_dim = 1;
        data_y_dim = 2;
        data_z_dim = 3;
        data_detector_dim = 4;
        data_X_dim = 5;
        data_Y_dim = 6;
    end
else
    if dimensions >= 4 %5D dataset
        data_detector_dim = -1; %No detector dimensions
        data_x_dim = 1;
        data_y_dim = 2;
        data_X_dim = 3;
        data_Y_dim = 4;
        data_z_dim = 5;
    else
        data_detector_dim = -1; %No detector dimensions
        data_x_dim = 1;
        data_y_dim = 2;
        data_z_dim = 3;
        data_X_dim = 4;
        data_Y_dim = 5;
    end
end
dims = [data_x_dim, data_y_dim, data_z_dim, data_detector_dim, data_X_dim, data_Y_dim]; %predefined data dimensions


% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
reciprocal_simulations = [21, 22, 31, 41, 51, 71, 81, 91];
%realspace_simulations = [11, 12, 32, 42, 52, 72, 82, 92];
if ismember(results_struct.input.simulation_type, reciprocal_simulations)
    units = '$Å^{-1}$';
else
    units = 'Å';
end
    
if dimensions >= data_x_dim
    axes_struct.x = struct();
    axes_struct.x.name = 'x';
    axes_struct.x.scale = results_struct.dx;
    axes_struct.x.offset = 0;
    axes_struct.x.units = units;
    axes_struct.x.size = int64(size(results_struct.images, data_x_dim));
    axes_struct.x.navigate = false;
end

    
if dimensions >= data_y_dim
    axes_struct.y = struct();
    axes_struct.y.name = 'y';
    axes_struct.y.scale = results_struct.dy;
    axes_struct.y.offset = 0;
    axes_struct.y.units = units;
    axes_struct.y.size = int64(size(results_struct.images, data_y_dim));
    axes_struct.y.navigate = false;
end

if dimensions >= data_X_dim
    axes_struct.X = struct();
    axes_struct.X.name = 'X';
    axes_struct.X.scale = results_struct.xs(2) - results_struct.xs(1);
    axes_struct.X.offset = results_struct.xs(1);
    axes_struct.X.units = 'Å';
    axes_struct.X.size = int64(size(results_struct.images, data_X_dim));
    axes_struct.X.navigate = true;
end

if dimensions >= data_Y_dim
    axes_struct.Y = struct();
    axes_struct.Y.name = 'Y';
    axes_struct.Y.scale = results_struct.ys(2) - results_struct.ys(1);
    axes_struct.Y.offset = results_struct.ys(1);
    axes_struct.Y.units = 'Å';
    axes_struct.Y.size = int64(size(results_struct.images, data_Y_dim));
    axes_struct.Y.navigate = true;
end

if dimensions > max(dims)
    for dim=5:dimensions
        axes_struct.(spintf('axes%i', dim)) = struct();
        axes_struct.(spintf('axes%i', dim)).name = sprintf('axes%i', dim);
        axes_struct.(spintf('axes%i', dim)).scale = 1;
        axes_struct.(spintf('axes%i', dim)).offset = 0;
        axes_struct.(spintf('axes%i', dim)).units = '';
        axes_struct.(spintf('axes%i', dim)).size = int64(size(results_struct.images, dim));
        axes_struct.(spintf('axes%i', dim)).navigate = true;
    end
end

if dimensions >= data_z_dim
    axes_struct.z = struct();
    axes_struct.z.name = 'z';
    axes_struct.z.scale = results_struct.thick(2) - results_struct.thick(1);%results_struct.input.spec_dz;
    axes_struct.z.offset = results_struct.thick(1);
    axes_struct.z.units = 'Å';
    axes_struct.z.size = int64(size(results_struct.images, data_z_dim));
    axes_struct.z.navigate = true;
    
end

if dimensions >= data_detector_dim && data_detector_dim > 0
    axes_struct.detectors = struct();
    axes_struct.detectors.name = 'detectors';
    axes_struct.detectors.scale = 1;
    axes_struct.detectors.offset = 0;
    axes_struct.detectors.units = '';
    axes_struct.detectors.size = int64(size(results_struct.images, data_detector_dim));
    axes_struct.detectors.navigate = true;
end

end

