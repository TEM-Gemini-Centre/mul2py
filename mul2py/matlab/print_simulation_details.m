function [details] = print_simulation_details(input_multislice, varargin)
    details = "";
    default_MULTEM_path = '/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM';
    default_print_parser = 0;
    
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && x>=0;
    validStrChar = @(x) ischar(x) || isstring(x);
    
    p = inputParser;
    addRequired(p, "input_multislice", @isstruct);
    addParameter(p, "print_parser", default_print_parser, validScalarPosNum);
    addParameter(p, "MULTEM_path", default_MULTEM_path, validStrChar);
    parse(p, input_multislice, varargin{:});
    
    if p.Results.print_parser
        fprintf("Parser:\n")
        disp(p)

        fprintf("Parser results:\n")
        disp(p.Results)
    end
    
    addpath(char(sprintf("%s/crystalline_materials", p.Results.MULTEM_path)));
    addpath(char(sprintf("%s/matlab_functions", p.Results.MULTEM_path)));
    addpath(char(sprintf("%s/mex_bin", p.Results.MULTEM_path)));
    
    %% Specimen:
    details = details + sprintf("Specimen dimensions:\n\tLx: %d Å\n\tLy: %d Å\n\tLz: %d Å\n", input_multislice.spec_lx, input_multislice.spec_ly, input_multislice.spec_lz);
    
    %% Slicing:
    details = details + sprintf("Slice thickness: %d Å\nSlices: %i\nStart: %d Å\nStop: %d Å\n", input_multislice.spec_dz, length(input_multislice.thick), input_multislice.thick(1), input_multislice.thick(end));
    
    %% Wavelength
    emass = 510.99906;		% electron rest mass in keV
    hc = 12.3984244;		% Planck's const x speed of light	
    lambda = hc/sqrt(input_multislice.E_0*(2.0*emass + input_multislice.E_0)); %wavelength of electron
    details = details + sprintf("Wavelength: %d Å\n", lambda);
    
    %% Bandwidth
    if input_multislice.bwl
        bandwidth = 2/3;
    else
        bandwidth = 1;
    end
    details = details + sprintf("Bandwidth: %d\n", bandwidth);
    
    %% Potential resolution
    dx = input_multislice.spec_lx / input_multislice.nx;
    dy = input_multislice.spec_ly / input_multislice.ny;
    details = details +sprintf("Potential resolution (%i, %i):\n\tx: %d Å \n\ty: %d Å\n", input_multislice.nx, input_multislice.ny, dx, dy);
    
    %% Reciprocal resolution
    dkx = 1 / input_multislice.spec_lx;
    dky = 1 / input_multislice.spec_ly;
    dax = rAng_2_mrad(dkx, input_multislice.E_0);
    day = rAng_2_mrad(dky, input_multislice.E_0);
    details = details +sprintf("Reciprocal resolution:\n\tx: %d mrad (%d 1/Å) \n\ty: %d mrad (%d 1/Å)\n", dax, dkx, day, dky);
    
    %% Maximum scattering angles
    %Maximum scattering vectors
    kx_max = input_multislice.nx * dkx * bandwidth / 2; %Divide by 2 because half-angles
    ky_max = input_multislice.ny * dky * bandwidth / 2;
    
    max_scattering_angle_x = rAng_2_mrad(kx_max, input_multislice.E_0);
    max_scattering_angle_y = rAng_2_mrad(ky_max, input_multislice.E_0);
    details = details + sprintf("Maximum scattering angles:\n\tx: %d mrad (%d 1/Å)\n\ty: %d mrad (%d 1/Å)\n", max_scattering_angle_x, kx_max, max_scattering_angle_y, ky_max);    
    
    