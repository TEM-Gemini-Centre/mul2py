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

function [input_multislice] = JEM2100F_EWRS_setup(model_path, alpha, varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%% Argument Parsing %%%%%%%%%%%%%%%%%%%%%%%%
    default_mode = 'converged';
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
    default_print_parser = 0;
    default_MULTEM_path = '/lustre1/projects/itea_lille-nv-fys-tem/MULTEM/MULTEM';
    
    p = inputParser;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    addRequired(p, 'model_path', @isstring);
    addRequired(p, 'alpha', validScalarPosNum);
    addParameter(p, 'mode', default_mode, @isstring);
    addParameter(p, 'nx', default_nx, validScalarPosNum);
    addParameter(p, 'ny', default_ny, validScalarPosNum);
    addParameter(p, 'bwl', default_bwl, validScalarPosNum);
    addParameter(p, 'x', default_x, validScalarPosNum);
    addParameter(p, 'y', default_y, validScalarPosNum);
    addParameter(p, 'E0', default_E0, validScalarPosNum);
    addParameter(p, 'phonons', default_phonons, validScalarPosNum);
    addParameter(p, 'thick_type', default_thick_type, validScalarPosNum);
    addParameter(p, 'thicknesses', default_thicknesses);
    addParameter(p, 'defocus', default_defocus);
    addParameter(p, 'print_parser', default_print_parser);
    addParameter(p, 'MULTEM_path', default_MULTEM_path, @isstring);
    parse(p, model_path, alpha, varargin{:});
    
    if p.Results.print_parser
        fprintf('Parser:\n')
        disp(p)

        fprintf('Parser results:\n')
        disp(p.Results)
    end
    
    addpath(sprintf('%s/crystalline_materials', p.Results.MULTEM_path));
    addpath(sprintf('%s/matlab_functions', p.Results.MULTEM_path));
    addpath(sprintf('%s/mex_bin', p.Results.MULTEM_path));

    %%%%%%%%%%%%%%%%%% Load multem default parameter %%%%%%%%$$%%%%%%%%%
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
        fprintf('Error when setting specimen crystal parameters "na", "nb", and "nc":\n\t%s\nSetting values based on specimen unit cell and simulation cell size. \nNB! only makes sense for cubic crystals with a, b, and c oriented along x, y, and z, respectively.\n', ME.message)
        input_multislice.spec_cryst_na = input_multislice.spec_lx / input_multislice.spec_cryst_a;
        input_multislice.spec_cryst_nb = input_multislice.spec_ly / input_multislice.spec_cryst_b;
        input_multislice.spec_cryst_nc = input_multislice.spec_lz / input_multislice.spec_cryst_c;
    end
    
    %%%%%%%%%%%%%%%%%%%%%% Adjust beam x-y ? %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isnan(p.Results.x)
    	x = input_multislice.spec_lx / 2;
    else
        x = p.Results.x;
    end
    if isnan(p.Results.y)
    	y = input_multislice.spec_ly / 2;
    else
    	y = p.Results.y;
    end

    %%%%%%%%%%%%%%%%%%%%%% Specimen thickness %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.thick_type = p.Results.thick_type;                     % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
    if p.Results.thick_type == 2
        if length(p.Results.thicknesses) == 1
            if p.Results.thicknesses == 0 %Use the slice thicknesses
                input_multislice.thick = (0:input_multislice.spec_dz:input_multislice.spec_lz+input_multislice.spec_dz);
            else %use the provided thickness
                input_multislice.thick = p.Results.thicknesses;
            end
        else %use the provided thicknesses
            input_multislice = p.Results.thicknesses;
        end
    end

    %%%%%%%%%%%%%%%%%%%% Set simulation experiment %%%%%%%%%%%%%%%%%%%%%
    % eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
    % eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
    input_multislice.simulation_type = 52;

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
    input_multislice.ny = p.Results.nx;
    input_multislice.bwl = p.Results.bwl;                            % Band-width limit, 1: true, 0:false

    %%%%%%%%%%%%%%%%%%%% Microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.E_0 = p.Results.E0;                          % Acceleration Voltage (keV)
    input_multislice.theta = 0.0;                        % Till ilumination (�)
    input_multislice.phi = 0.0;                          % Till ilumination (�)

    %%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.illumination_model = 1;             % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
    input_multislice.temporal_spatial_incoh = 1;         % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(p.Results.mode, 'converged')
        input_multislice.iw_type = 2;                        % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
    elseif strcmp(p.Results.mode, 'plane')
        input_multislice.iw_type = 1;
    else
        input_multislice.iw_type = 4;
    end
    input_multislice.iw_psi = 0;    % user define incident wave
    input_multislice.iw_x = x;  % x position 
    input_multislice.iw_y = y;  % y position

    %%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.cond_lens_m = 0;                   % Vortex momentum
    input_multislice.cond_lens_c_10 = p.Results.defocus;             % Defocus (�)
    input_multislice.cond_lens_c_30 = 1.0;              % Third order spherical aberration (mm)
    input_multislice.cond_lens_c_50 = 0.00;             % Fifth order spherical aberration (mm)
    input_multislice.cond_lens_c_12 = 0.0;              % Twofold astigmatism (�)
    input_multislice.cond_lens_phi_12 = 0.0;            % Azimuthal angle of the twofold astigmatism (�)
    input_multislice.cond_lens_c_23 = 0.0;              % Threefold astigmatism (�)
    input_multislice.cond_lens_phi_23 = 0.0;            % Azimuthal angle of the threefold astigmatism (�)
    input_multislice.cond_lens_inner_aper_ang = 0.0;    % Inner aperture (mrad) 
    input_multislice.cond_lens_outer_aper_ang = p.Results.alpha;   % Outer aperture (mrad)

    %%%%%%%%% defocus spread function %%%%%%%%%%%%
    dsf_sigma = il_iehwgd_2_sigma(32); % from defocus spread to standard deviation
    input_multislice.cond_lens_dsf_sigma = dsf_sigma;   % standard deviation (�)
    input_multislice.cond_lens_dsf_npoints = 5;         % # of integration points. It will be only used if illumination_model=4

    %%%%%%%%%% source spread function %%%%%%%%%%%%
    ssf_sigma = il_hwhm_2_sigma(0.45); % half width at half maximum to standard deviation
    input_multislice.cond_lens_ssf_sigma = ssf_sigma;  	% standard deviation: For parallel ilumination(�^-1); otherwise (�)
    input_multislice.cond_lens_ssf_npoints = 4;         % # of integration points. It will be only used if illumination_model=4

    %%%%%%%%% zero defocus reference %%%%%%%%%%%%
    input_multislice.cond_lens_zero_defocus_type = 1;   % eZDT_First = 1, eZDT_User_Define = 2
    input_multislice.cond_lens_zero_defocus_plane = 0;