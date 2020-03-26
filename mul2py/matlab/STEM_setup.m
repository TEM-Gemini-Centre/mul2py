% Function for setting up MULTEM CBED simulations for the JEM2100F 
% microscope at the TEM Gemini centre.
% Based on the CBED example written by Ivan Lobato:
%
% output_multislice = il_MULTEM(system_conf, input_multislice) perform TEM simulation
% 
% Exit wave real space (EWRS) simulation
% 
% All parameters of the input_multislice structure are explained in multem_default_values()
% 
% Copyright 2017 Ivan Lobato <Ivanlh20@gmail.com>

function [input_multislice] = STEM_setup(model_path, alpha, collection_angles, varargin)
    
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
    default_MULTEM_path = '/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM';
    
    p = inputParser;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
    validScalarNum = @(x) isnumeric(x) && isscalar(x);
    validStrChar = @(x) ischar(x) || isstring(x);
    
    addRequired(p, "model_path", validScalarPosNum);
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
    clear input_multislice
    input_multislice = multem_default_values();          % Load default values;

    %%%%%%%%%%%%%%%%%%%%%%% Specimen information %%%%%%%%%%%%%%%%%%%%%%%
    load(p.Results.model_path);
    %Update `input_multislice`
    input_multislice.spec_atoms = spec_atoms;
    input_multislice.spec_lx = spec_lx;
    input_multislice.spec_ly = spec_ly;
    input_multislice.spec_lz = spec_lz;
    input_multislice.spec_dz = spec_dz;
    input_multislice.spec_cryst_a = a;
    input_multislice.spec_cryst_b = b;
    input_multislice.spec_cryst_c = c;
    try
        input_multislice.spec_cryst_na = na;
        input_multislice.spec_cryst_nb = nb;
        input_multislice.spec_cryst_nc = nc;
    catch ME
        fprintf("Error when setting specimen crystal parameters 'na', 'nb', and 'nc':\n\t%s\nSetting values based on specimen unit cell and simulation cell size. \nNB! only makes sense for cubic crystals with a, b, and c oriented along x, y, and z, respectively.\n", ME.message)
        input_multislice.spec_cryst_na = input_multislice.spec_lx / input_multislice.spec_cryst_a;
        input_multislice.spec_cryst_nb = input_multislice.spec_ly / input_multislice.spec_cryst_b;
        input_multislice.spec_cryst_nc = input_multislice.spec_lz / input_multislice.spec_cryst_c;
    end
        

    %%%%%%%%%%%%%%%%%%%%%% Specimen thickness %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.thick_type = p.Results.thick_type;                     % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
    if p.Results.thick_type == 2
        if length(p.Results.thicknesses) == 1
            if p.Results.thicknesses == 0 %Use the slice thicknesses
                input_multislice.thick = (0:input_multislice.spec_dz:input_multislice.spec_lz-input_multislice.spec_dz);
            else %use the provided thickness
                input_multislice.thick = p.Results.thicknesses;
            end
        else %use the provided thicknesses
            input_multislice = p.Results.thicknesses;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%% Scanning region %%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.scanning_type = 2; % eST_Line = 1, eST_Area = 2
    input_multislice.scanning_periodic = 0;     % 1: true, 0:false (periodic boundary conditions)
    input_multislice.scanning_ns = p.Results.scanning_ns;
    
    if isnan(p.Results.center_x)
    	center_x = input_multislice.spec_lx / 2;
    else
        center_x = p.Results.center_x;
    end
    if isnan(p.Results.center_y)
    	center_y = input_multislice.spec_ly / 2;
    else
    	center_y = p.Results.center_y;
    end
    
    if isnan(p.Results.scan_width)
        scan_width = input_multislice.spec_cryst_a;
    else
        scan_width = p.Results.scan_width;
    end
    
    if isnan(p.Results.scan_height)
        scan_height = input_multislice.spec_cryst_b;
    else
        scan_height = p.Results.scan_height;
    end
    
    input_multislice.scanning_x0 = center_x - scan_width/2; 
    input_multislice.scanning_y0 = center_y - scan_height/2;
    input_multislice.scanning_xe = input_multislice.scanning_x0 + scan_width;
    input_multislice.scanning_ye = input_multislice.scanning_y0 + scan_height;

    %%%%%%%%%%%%%%%%%%%% Set simulation experiment %%%%%%%%%%%%%%%%%%%%%
    % eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
    % eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
    input_multislice.simulation_type = 11;

    %%%%%%%%%%%%%% Electron-Specimen interaction model %%%%%%%%%%%%%%%%%
    input_multislice.interaction_model = 1;              % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
    input_multislice.potential_type = 6;                 % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

    %%%%%%%%%%%%%%%%%%%%%%% Potential slicing %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.potential_slicing = 2;              % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4

    %%%%%%%%%%%%%%% Electron-Phonon interaction model %%%%%%%%%%%%%%%%%%
    input_multislice.pn_model = 3;                       % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
    input_multislice.pn_coh_contrib = 0;
    input_multislice.pn_single_conf = 0;                 % 1: true, 0:false (extract single configuration)
    input_multislice.pn_nconf = p.Results.phonons;                      % true: specific phonon configuration, false: number of frozen phonon configurations
    input_multislice.pn_dim = 110;                       % phonon dimensions (xyz)
    input_multislice.pn_seed = 300183;                   % Random seed(frozen phonon)

    %%%%%%%%%%%%%%%%%%%%%% x-y sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.nx = p.Results.nx;
    input_multislice.ny = p.Results.ny;
    input_multislice.bwl = p.Results.bwl;                            % Band-width limit, 1: true, 0:false

    %%%%%%%%%%%%%%%%%%%% Microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.E_0 = p.Results.E0;                          % Acceleration Voltage (keV)
    input_multislice.theta = 0.0;                        % Till ilumination (ï¿½)
    input_multislice.phi = 0.0;                          % Till ilumination (ï¿½)

    %%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.illumination_model = 1;             % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
    input_multislice.temporal_spatial_incoh = 1;         % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

    %%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.cond_lens_m = 0;                   % Vortex momentum
    input_multislice.cond_lens_inner_aper_ang = 0.0;    % Inner aperture (mrad) 
    input_multislice.cond_lens_outer_aper_ang = p.Results.alpha;   % Outer aperture (mrad)
    
    % Set aberrations
    if strcmp(p.Results.instrument, "")
        aberrations = nan;
    elseif strcmp(p.Results.instrument, "ARM200F")
        aberrations = ARM200F_aberrations();
    elseif strcmp(p.Results.instrument, "2100F")
        aberrations = JEM2100F_aberrations();
    else
        fprintf("Could not understand instrument %s, using default aberrations", p.Results.instrument);
        aberrations = nan;
    end
    
    if isstruct(aberrations)
        input_multislice.cond_lens_c_12 = aberrations.cond_lens_c_12;                         % [A1]      2-fold astigmatism (Å)
        input_multislice.cond_lens_c_phi_12 = aberrations.cond_lens_phi_12;                   % [phi_A1]	Azimuthal angle of 2-fold astigmatism (deg)

        input_multislice.cond_lens_c_21 = aberrations.cond_lens_c_21;                         % [B2]      Axial coma (Å)
        input_multislice.cond_lens_c_phi_21 = aberrations.cond_lens_phi_21;                   % [phi_B2]	Azimuthal angle of axial coma (deg)

        input_multislice.cond_lens_c_23 = aberrations.cond_lens_c_23;                         % [A2]      3-fold astigmatism (Å)
        input_multislice.cond_lens_c_phi_23 = aberrations.cond_lens_phi_23;                   % [phi_A2]	Azimuthal angle of 3-fold astigmatism (deg)

        input_multislice.cond_lens_c_30 = aberrations.cond_lens_c_30;                         % [C3] 		3rd order spherical aberration (mm)

        input_multislice.cond_lens_c_32 = aberrations.cond_lens_c_32;                         % [S3]      Axial star aberration (Å)
        input_multislice.cond_lens_c_phi_32 = aberrations.cond_lens_phi_32;                   % [phi_S3]	Azimuthal angle of axial star aberration (deg)

        input_multislice.cond_lens_c_34 = aberrations.cond_lens_c_34;                         % [A3]      4-fold astigmatism (Å)
        input_multislice.cond_lens_c_phi_34 = aberrations.cond_lens_phi_34;                   % [phi_A3]	Azimuthal angle of 4-fold astigmatism (deg)

        input_multislice.cond_lens_c_41 = aberrations.cond_lens_c_41;                         % [B4]      4th order axial coma (Å)
        input_multislice.cond_lens_c_phi_41 = aberrations.cond_lens_phi_41;                   % [phi_B4]	Azimuthal angle of 4th order axial coma (deg)

        input_multislice.cond_lens_c_43 = aberrations.cond_lens_c_43;                         % [D4]      3-lobe aberration (Å)
        input_multislice.cond_lens_c_phi_43 = aberrations.cond_lens_phi_43;                   % [phi_D4]	Azimuthal angle of 3-lobe aberration (deg)

        input_multislice.cond_lens_c_45 = aberrations.cond_lens_c_45;                         % [A4]      5-fold astigmatism (Å)
        input_multislice.cond_lens_c_phi_45 = aberrations.cond_lens_phi_45;                   % [phi_A4]	Azimuthal angle of 5-fold astigmatism (deg)

        input_multislice.cond_lens_c_50 = aberrations.cond_lens_c_50;                         % [C5]      5th order spherical aberration (mm)
        input_multislice.cond_lens_c_52 = aberrations.cond_lens_c_52;                         % [S5]      5th order axial star aberration (?)
        input_multislice.cond_lens_c_phi_52 = aberrations.cond_lens_phi_52;                   % [phi_S5]	Azimuthal angle of 5th order axial star aberration (?)
        input_multislice.cond_lens_c_54 = aberrations.cond_lens_c_54;                         % [R5]      5th order rosette aberration (?)
        input_multislice.cond_lens_c_phi_54 = aberrations.cond_lens_phi_54;                   % [phi_R5]	Azimuthal angle of 5th order rosette aberration (?)
        input_multislice.cond_lens_c_56 = aberrations.cond_lens_c_56;                         % [A5]      6-fold astigmatism (?)
        input_multislice.cond_lens_c_phi_56 = aberrations.cond_lens_phi_56;                   % [phi_A5]	Azimuthal angle of 6-fold astigmatism (?)
    end
    
    if isnan(p.Results.defocus)
    	defocus = il_scherzer_defocus(input_multislice.E_0, input_multislice.obj_lens_c_30);
    else
        defocus = p.Results.defocus;
    end
    input_multislice.obj_lens_c_10 = defocus;

    %%%%%%%%% defocus spread function %%%%%%%%%%%%
    dsf_sigma = il_iehwgd_2_sigma(32); % from defocus spread to standard deviation
    input_multislice.cond_lens_dsf_sigma = dsf_sigma;   % standard deviation (ï¿½)
    input_multislice.cond_lens_dsf_npoints = 5;         % # of integration points. It will be only used if illumination_model=4

    %%%%%%%%%% source spread function %%%%%%%%%%%%
    ssf_sigma = il_hwhm_2_sigma(0.45); % half width at half maximum to standard deviation
    input_multislice.cond_lens_ssf_sigma = ssf_sigma;  	% standard deviation: For parallel ilumination(ï¿½^-1); otherwise (ï¿½)
    input_multislice.cond_lens_ssf_npoints = 4;         % # of integration points. It will be only used if illumination_model=4

    %%%%%%%%% zero defocus reference %%%%%%%%%%%%
    input_multislice.cond_lens_zero_defocus_type = 1;   % eZDT_First = 1, eZDT_User_Define = 2
    input_multislice.cond_lens_zero_defocus_plane = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Detector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.detector.type = 1;  % eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3
    %Set collection angles (in mrad) of circular detectors
    for d=1:length(collection_angles)
        input_multislice.detector.cir(d) = collection_angles(d);
    end
    
    if p.Results.print_details
        fprintf("**** Set up MULTEM STEM simulation for instrument '%s' ****\n\n%s\n", p.Results.instrument, print_simulation_details(input_multislice, "MULTEM_path", p.Results.MULTEM_path))
    end