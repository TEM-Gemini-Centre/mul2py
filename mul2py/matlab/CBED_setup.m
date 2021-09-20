% Function for setting up MULTEM SCBED simulations for the JEM2100F
% microscope at the TEM Gemini centre.
% Based on the SCBED example written by Ivan Lobato:
%
% output_multislice = il_MULTEM(system_conf, input_multislice) perform TEM simulation
% 
% 
% All parameters of the input_multislice structure are explained in multem_default_values()
% 
% Copyright 2017 Ivan Lobato <Ivanlh20@gmail.com>

function [input_multem] = CBED_setup(model_path, alpha, varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%% Argument Parsing %%%%%%%%%%%%%%%%%%%%%%%%
    default_defocus = 0;
    default_nx = 1024;
    default_ny = 1024;
    default_bwl = 1;
    default_x = nan;
    default_y = nan;
    default_E0 = 200;
    default_phonons = 20;
    default_thick_type = 2;
    default_thicknesses = 0;
    default_instrument = "";
    default_print_parser = 0;
    default_print_details = 1;
    default_MULTEM_path = "/cluster/projects/itea_lille-nv-fys-tem/repositories/multem";
    
    p = inputParser;
    p.KeepUnmatched = true;
    
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
    validScalarNum = @(x) isnumeric(x) && isscalar(x);
    validStrChar = @(x) ischar(x) || isstring(x);
    
    addRequired(p, "model_path", validStrChar);
    addRequired(p, "alpha", validScalarPosNum);
    addParameter(p, "nx", default_nx, validScalarPosNum);
    addParameter(p, "ny", default_ny, validScalarPosNum);
    addParameter(p, "bwl", default_bwl, validScalarPosNum);
    addParameter(p, "x", default_x, validScalarNum);
    addParameter(p, "y", default_y, validScalarNum);
    addParameter(p, "E0", default_E0, validScalarPosNum);
    addParameter(p, "phonons", default_phonons, validScalarPosNum);
    addParameter(p, "thick_type", default_thick_type, validScalarPosNum);
    addParameter(p, "thicknesses", default_thicknesses);
    addParameter(p, "defocus", default_defocus, validScalarNum);
    addParameter(p, "instrument", default_instrument, validStrChar);
    addParameter(p, "print_parser", default_print_parser, validScalarPosNum);
    addParameter(p, "print_details", default_print_details, validScalarPosNum);
    addParameter(p, "MULTEM_path", default_MULTEM_path, validStrChar);
    parse(p, model_path, alpha, varargin{:});
    
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

    %%%%%%%%%%%%%%%%%%%% Set simulation experiment %%%%%%%%%%%%%%%%%%%%%
    % eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
    % eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
    input_multem.simulation_type = 21;

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
    input_multem.illumination_model = 4;%1;             % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
    input_multem.temporal_spatial_incoh = 1;         % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.iw_type = 4;                        % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
    input_multem.iw_psi = read_psi_0_multem(input_multem.nx, input_multem.ny);    % user define incident wave
    
    %%%%%%%%%%%%%%%%%%%%%% Adjust beam x-y  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isnan(p.Results.x)
    	x = input_multem.spec_lx / 2;
    else
        x = p.Results.x;
    end
    if isnan(p.Results.y)
    	y = input_multem.spec_ly / 2;
    else
    	y = p.Results.y;
    end
    input_multem.iw_x = x;  % x position 
    input_multem.iw_y = y;  % y position

    %%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.cond_lens_m = 0;                   % Vortex momentum
    
    input_multem.cond_lens_inner_aper_ang = 0.0;    % Inner aperture (mrad) 
    input_multem.cond_lens_outer_aper_ang = p.Results.alpha;   % Outer aperture (mrad)
    
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
        input_multem = set_aberrations(input_multem, aberrations);
    end
    
    defocus = p.Results.defocus;
    input_multem.cond_lens_c_10 = defocus;

    %%%%%%%%% defocus spread function %%%%%%%%%%%%
    dsf_sigma = ilc_iehwgd_2_sigma(32); % from defocus spread to standard deviation
    input_multem.cond_lens_ti_a = 1.0;                          % Height proportion of a normalized Gaussian [0, 1]
    input_multem.cond_lens_ti_sigma = dsf_sigma;                % Standard deviation of the defocus spread for the Gaussian component (�)
    input_multem.cond_lens_ti_beta = 0.0;                 		% Standard deviation of the defocus spread for the Exponential component (�)
    input_multem.cond_lens_ti_npts = 4;                         % Number of integration points. It will be only used if illumination_model=4

    %%%%%%%%%% source spread function %%%%%%%%%%%%
    ssf_sigma = ilc_hwhm_2_sigma(0.45); % half width at half maximum to standard deviation
    input_multem.cond_lens_si_a = 1.0;                          % Height proportion of a normalized Gaussian [0, 1]
    input_multem.cond_lens_si_sigma = ssf_sigma;                % Standard deviation of the source spread function for the Gaussian component: For parallel ilumination(�^-1); otherwise (�)
    input_multem.cond_lens_si_beta = 0.0;                 		% Standard deviation of the source spread function for the Exponential component: For parallel ilumination(�^-1); otherwise (�)
    input_multem.cond_lens_si_rad_npts = 4;                     % Number of radial integration points. It will be only used if illumination_model=4
    input_multem.cond_lens_si_azm_npts = 4;                     % Number of radial integration points. It will be only used if illumination_model=4

    %%%%%%%%% zero defocus reference %%%%%%%%%%%%
    input_multem.cond_lens_zero_defocus_type = 1;   % eZDT_First = 1, eZDT_User_Define = 2
    input_multem.cond_lens_zero_defocus_plane = 0;
    
    if p.Results.print_details
        fprintf("**** Set up MULTEM CBED simulation for instrument '%s' ****\n\n%s\n", p.Results.instrument, print_simulation_details(input_multem, "MULTEM_path", p.Results.MULTEM_path))
    end