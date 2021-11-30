% Function for setting up MULTEM SCBED simulations for the JEM2100F
% microscope at the TEM Gemini centre.
% Based on the SCBED example written by Ivan Lobato:
%
% output_multislice = il_MULTEM(system_conf, input_multislice) perform TEM simulation
% 
% Exit wave real space (EWRS) simulation
% 
% All parameters of the input_multislice structure are explained in multem_default_values()
% 
% Copyright 2017 Ivan Lobato <Ivanlh20@gmail.com>

function [input_multem] = STEM_setup(model_path, alpha, collection_angles, varargin)
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Argument Parsing %%%%%%%%%%%%%%%%%%%%%%%%
    default_defocus = 0;
    default_nx = 1024;
    default_ny = 1024;
    default_bwl = 1;
    default_center_x = nan;
    default_center_y = nan;
    default_scanning_ns = 10;
    default_scan_width = nan;
    default_scan_height = nan;
    default_E0 = 200;
    default_phonons = 20;
    default_thick_type = 2;
    default_thicknesses = 0;
    default_instrument = "";
    default_print_parser = 0;
    default_print_details = 1;
    default_scanning_periodic = 1;
    default_MULTEM_path = '/cluster/projects/itea_lille-nv-fys-tem/repositories/multem';
    
    p = inputParser;
    p.KeepUnmatched = true;

    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
    validScalarNum = @(x) isnumeric(x) && isscalar(x);
    validStrChar = @(x) ischar(x) || isstring(x);
    
    addRequired(p, "model_path", validStrChar);
    addRequired(p, "alpha", validScalarPosNum);
    addRequired(p, "collection_angles", @isstruct)
    addParameter(p, "nx", default_nx, validScalarPosNum);
    addParameter(p, "ny", default_ny, validScalarPosNum);
    addParameter(p, "bwl", default_bwl, validScalarPosNum);
    addParameter(p, "center_x", default_center_x, validScalarNum);
    addParameter(p, "center_y", default_center_y, validScalarNum);
    addParameter(p, "scan_width", default_scan_width, validScalarPosNum);
    addParameter(p, "scan_height", default_scan_height, validScalarPosNum);
    addParameter(p, "scanning_ns", default_scanning_ns, validScalarPosNum);
    addParameter(p, "E0", default_E0, validScalarPosNum);
    addParameter(p, "phonons", default_phonons, validScalarPosNum);
    addParameter(p, "thick_type", default_thick_type, validScalarPosNum);
    addParameter(p, "thicknesses", default_thicknesses);
    addParameter(p, "defocus", default_defocus, validScalarNum);
    addParameter(p, "instrument", default_instrument, validStrChar);
    addParameter(p, "print_parser", default_print_parser, validScalarPosNum);
    addParameter(p, "print_details", default_print_details, validScalarPosNum);
    addParameter(p, "MULTEM_path", default_MULTEM_path, validStrChar);
    addParameter(p, "scanning_periodic", default_scanning_periodic, validScalarPosNum);
    parse(p, model_path, alpha, collection_angles, varargin{:});
    
    if p.Results.print_parser
        fprintf("Parser:\n")
        disp(p)
        
        fprintf("Parser results:\n")
        disp(p.Results)
    end
    
    addpath(char(sprintf("%s/crystalline_materials", p.Results.MULTEM_path)));
    addpath(char(sprintf("%s/matlab_functions", p.Results.MULTEM_path)));
    addpath(char(sprintf("%s/mex_bin", p.Results.MULTEM_path)));
    
    
    %%%%%%%%%%%%%%%%%% Load multem default parameter %%%%%%%%$$%%%%%%%%%
    clear input_multem
    input_multem = multem_input.parameters;          % Load default values;

    %%%%%%%%%%%%%%%%%%%%%%% Specimen information %%%%%%%%%%%%%%%%%%%%%%%
    load(p.Results.model_path);
    %Update `input_multem`
    input_multem.spec_atoms = spec_atoms;
    input_multem.spec_lx = spec_lx;
    input_multem.spec_ly = spec_ly;
    input_multem.spec_lz = spec_lz;
    input_multem.spec_dz = spec_dz;
    input_multem.spec_cryst_a = a;
    input_multem.spec_cryst_b = b;
    input_multem.spec_cryst_c = c;
    try
        input_multem.spec_cryst_na = na;
        input_multem.spec_cryst_nb = nb;
        input_multem.spec_cryst_nc = nc;
    catch ME
        fprintf("Error when setting specimen crystal parameters 'na', 'nb', and 'nc':\n\t%s\nSetting values based on specimen unit cell and simulation cell size. \nNB! only makes sense for cubic crystals with a, b, and c oriented along x, y, and z, respectively.\n", ME.message)
        input_multem.spec_cryst_na = input_multem.spec_lx / input_multem.spec_cryst_a;
        input_multem.spec_cryst_nb = input_multem.spec_ly / input_multem.spec_cryst_b;
        input_multem.spec_cryst_nc = input_multem.spec_lz / input_multem.spec_cryst_c;
    end
        

    %%%%%%%%%%%%%%%%%%%%%% Specimen thickness %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.thick_type = p.Results.thick_type;                     % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
    if p.Results.thick_type == 2
        if length(p.Results.thicknesses) == 1
            if p.Results.thicknesses == 0 %Use the slice thicknesses
                input_multem.thick = (0:input_multem.spec_dz:input_multem.spec_lz-input_multem.spec_dz);
            else %use the provided thickness
                input_multem.thick = p.Results.thicknesses;
            end
        else %use the provided thicknesses
            input_multem = p.Results.thicknesses;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%% Scanning region %%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.scanning_type = 2; % eST_Line = 1, eST_Area = 2
    input_multem.scanning_periodic = p.Results.scanning_periodic;     % 1: true, 0:false (periodic boundary conditions)
    input_multem.scanning_ns = p.Results.scanning_ns;
    
    if isnan(p.Results.center_x)
    	center_x = input_multem.spec_lx / 2;
    else
        center_x = p.Results.center_x;
    end
    if isnan(p.Results.center_y)
    	center_y = input_multem.spec_ly / 2;
    else
    	center_y = p.Results.center_y;
    end
    
    if isnan(p.Results.scan_width)
        scan_width = input_multem.spec_cryst_a;
    else
        scan_width = p.Results.scan_width;
    end
    
    if isnan(p.Results.scan_height)
        scan_height = input_multem.spec_cryst_b;
    else
        scan_height = p.Results.scan_height;
    end
    
    input_multem.scanning_x0 = center_x - scan_width/2; 
    input_multem.scanning_y0 = center_y - scan_height/2;
    input_multem.scanning_xe = input_multem.scanning_x0 + scan_width;
    input_multem.scanning_ye = input_multem.scanning_y0 + scan_height;

    %%%%%%%%%%%%%%%%%%%% Set simulation experiment %%%%%%%%%%%%%%%%%%%%%
    % eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
    % eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
    input_multem.simulation_type = 11;

    %%%%%%%%%%%%%% Electron-Specimen interaction model %%%%%%%%%%%%%%%%%
    input_multem.interaction_model = 1;              % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
    input_multem.potential_type = 6;                 % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

    %%%%%%%%%%%%%%%%%%%%%%% Potential slicing %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.potential_slicing = 2;              % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4

    %%%%%%%%%%%%%%% Electron-Phonon interaction model %%%%%%%%%%%%%%%%%%
    input_multem.pn_model = 3;                       % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
    input_multem.pn_coh_contrib = 0;
    input_multem.pn_single_conf = 0;                 % 1: true, 0:false (extract single configuration)
    input_multem.pn_nconf = p.Results.phonons;                      % true: specific phonon configuration, false: number of frozen phonon configurations
    input_multem.pn_dim = 110;                       % phonon dimensions (xyz)
    input_multem.pn_seed = 300183;                   % Random seed(frozen phonon)

    %%%%%%%%%%%%%%%%%%%%%% x-y sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.nx = p.Results.nx;
    input_multem.ny = p.Results.ny;
    input_multem.bwl = p.Results.bwl;                            % Band-width limit, 1: true, 0:false

    %%%%%%%%%%%%%%%%%%%% Microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.E_0 = p.Results.E0;                          % Acceleration Voltage (keV)
    input_multem.theta = 0.0;                        % Till ilumination (�)
    input_multem.phi = 0.0;                          % Till ilumination (�)

    %%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.illumination_model = 1;             % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
    input_multem.temporal_spatial_incoh = 1;         % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

    %%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.cond_lens_m = 0;                   % Vortex momentum
    input_multem.cond_lens_inner_aper_ang = 0.0;    % Inner aperture (mrad) 
    input_multem.cond_lens_outer_aper_ang = p.Results.alpha;   % Outer aperture (mrad)
    
    % Set aberrations
    if strcmp(p.Results.instrument, "")
        aberrations = no_aberrations();
    elseif strcmp(p.Results.instrument, "ARM200F")
        aberrations = ARM200F_aberrations();
    elseif strcmp(p.Results.instrument, "2100F")
        aberrations = JEM2100F_aberrations();
    else
        fprintf("Could not understand instrument %s, using default aberrations", p.Results.instrument);
        aberrations = no_aberrations();
    end
    
    if isstruct(aberrations)
        input_multem = set_aberrations(input_multem, aberrations);
    end
    
    defocus = p.Results.defocus;
    input_multem.cond_lens_c_10 = defocus;

    %%%%%%%%% defocus spread function %%%%%%%%%%%%
    dsf_sigma = ilc_iehwgd_2_sigma(32); % from defocus spread to standard deviation
    input_multem.cond_lens_ti_sigma = dsf_sigma;   % standard deviation (�)
    input_multem.cond_lens_ti_npts = 5;         % # of integration points. It will be only used if illumination_model=4

    %%%%%%%%%% source spread function %%%%%%%%%%%%
    ssf_sigma = ilc_hwhm_2_sigma(0.45); % half width at half maximum to standard deviation
    input_multem.cond_lens_si_sigma = ssf_sigma;  	% standard deviation: For parallel ilumination(�^-1); otherwise (�)
    input_multem.cond_lens_si_rad_npts = 4;         % # of integration points. It will be only used if illumination_model=4

    %%%%%%%%% zero defocus reference %%%%%%%%%%%%
    input_multem.cond_lens_zero_defocus_type = 1;   % eZDT_First = 1, eZDT_User_Define = 2
    input_multem.cond_lens_zero_defocus_plane = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Detector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.detector.type = 1;  % eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3
    %Set collection angles (in mrad) of circular detectors
    for d=1:length(collection_angles)
        input_multem.detector.cir(d) = collection_angles(d);
    end
    
    if p.Results.print_details
        fprintf("**** Set up MULTEM STEM simulation for instrument '%s' ****\n\n%s\n", p.Results.instrument, print_simulation_details(input_multem, "MULTEM_path", p.Results.MULTEM_path))
    end