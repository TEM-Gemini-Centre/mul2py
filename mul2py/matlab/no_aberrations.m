%ARM200F aberrations based on alignment 16/1/2017 (Adrian, Emil, Per Erik)
function [aberrations] = ARM200F_aberrations()
    clear aberrations
    %% CL aberrations
    aberrations.cond_lens_c_10 = 0.00;                         % [C1]      Defocus (�)
    
    aberrations.cond_lens_c_12 = 0.00;                   % [A1]      2-fold astigmatism (�)
    aberrations.cond_lens_phi_12 = 0.00;                       % [phi_A1]	Azimuthal angle of 2-fold astigmatism (deg)

    aberrations.cond_lens_c_21 = 0.00;                  % [B2]      Axial coma (�)
    aberrations.cond_lens_phi_21 = 0.00;                       % [phi_B2]	Azimuthal angle of axial coma (deg)
    
    aberrations.cond_lens_c_23 = 0.00;                  % [A2]      3-fold astigmatism (�)
    aberrations.cond_lens_phi_23 = 0.00;                       % [phi_A2]	Azimuthal angle of 3-fold astigmatism (deg)

    aberrations.cond_lens_c_30 = 0.00;                  % [C3] 		3rd order spherical aberration (mm)
    
    aberrations.cond_lens_c_32 = 0.00;                   % [S3]      Axial star aberration (�)
    aberrations.cond_lens_phi_32 = 0.00;                       % [phi_S3]	Azimuthal angle of axial star aberration (deg)
    
    aberrations.cond_lens_c_34 = 0.00;                   % [A3]      4-fold astigmatism (�)
    aberrations.cond_lens_phi_34 = 0.0;                        % [phi_A3]	Azimuthal angle of 4-fold astigmatism (deg)

    aberrations.cond_lens_c_41 = 0.00;                  % [B4]      4th order axial coma (�)
    aberrations.cond_lens_phi_41 = 0.00;                       % [phi_B4]	Azimuthal angle of 4th order axial coma (deg)
    
    aberrations.cond_lens_c_43 = 0.00;                   % [D4]      3-lobe aberration (�)
    aberrations.cond_lens_phi_43 = 0.00;                       % [phi_D4]	Azimuthal angle of 3-lobe aberration (deg)
    
    aberrations.cond_lens_c_45 = 0.00;                  % [A4]      5-fold astigmatism (�)
    aberrations.cond_lens_phi_45 = 0.00;                       % [phi_A4]	Azimuthal angle of 5-fold astigmatism (deg)

    aberrations.cond_lens_c_50 = 0.00;                         % [C5]      5th order spherical aberration (mm)
    aberrations.cond_lens_c_52 = 0.00;                         % [S5]      5th order axial star aberration (?)
    aberrations.cond_lens_phi_52 = 0.00;                       % [phi_S5]	Azimuthal angle of 5th order axial star aberration (?)
    aberrations.cond_lens_c_54 = 0.00;                         % [R5]      5th order rosette aberration (?)
    aberrations.cond_lens_phi_54 = 0.00;                       % [phi_R5]	Azimuthal angle of 5th order rosette aberration (?)
    aberrations.cond_lens_c_56 = 0.00;                         % [A5]      6-fold astigmatism (?)
    aberrations.cond_lens_phi_56 = 0.00;                       % [phi_A5]	Azimuthal angle of 6-fold astigmatism (?)

    %% OL aberrations
    aberrations.obj_lens_c_10 = 0.00;                          % [C1]      Defocus (�)
    
    aberrations.obj_lens_c_12 = 0.00;                      % [A1]      2-fold astigmatism (�)
    aberrations.obj_lens_phi_12 = 0.00;                        % [phi_A1]	Azimuthal angle of 2-fold astigmatism (deg)

    aberrations.obj_lens_c_21 = 0.00;                     % [B2]      Axial coma (�)
    aberrations.obj_lens_phi_21 = 0.00;                        % [phi_B2]	Azimuthal angle of axial coma (deg)
    
    aberrations.obj_lens_c_23 = 0.00;                      % [A2]      3-fold astigmatism (�)
    aberrations.obj_lens_phi_23 = 0.00;                        % [phi_A2]	Azimuthal angle of 3-fold astigmatism (deg)

    aberrations.obj_lens_c_30 = 0.00;                   % [C3] 		3rd order spherical aberration (mm)
    
    aberrations.obj_lens_c_32 = 0.00;                   % [S3]      Axial star aberration (�)
    aberrations.obj_lens_phi_32 = 0.00;                        % [phi_S3]	Azimuthal angle of axial star aberration (deg)
    
    aberrations.obj_lens_c_34 = 0.00;                   % [A3]      4-fold astigmatism (�)
    aberrations.obj_lens_phi_34 = 0.0;                         % [phi_A3]	Azimuthal angle of 4-fold astigmatism (deg)

    aberrations.obj_lens_c_41 = 0.00;                  % [B4]      4th order axial coma (�)
    aberrations.obj_lens_phi_41 = 0.00;                        % [phi_B4]	Azimuthal angle of 4th order axial coma (deg)
    
    aberrations.obj_lens_c_43 = 0.00;                  % [D4]      3-lobe aberration (�)
    aberrations.obj_lens_phi_43 = 0.00;                        % [phi_D4]	Azimuthal angle of 3-lobe aberration (deg)
    
    aberrations.obj_lens_c_45 = 0.00;                  % [A4]      5-fold astigmatism (�)
    aberrations.obj_lens_phi_45 = 0.00;                        % [phi_A4]	Azimuthal angle of 5-fold astigmatism (deg)

    aberrations.obj_lens_c_50 = 0.00;                          % [C5]      5th order spherical aberration (mm)
    aberrations.obj_lens_c_52 = 0.00;                          % [S5]      5th order axial star aberration (�)
    aberrations.obj_lens_phi_52 = 0.00;                        % [phi_S5]	Azimuthal angle of 5th order axial star aberration (deg)
    aberrations.obj_lens_c_54 = 0.00;                          % [R5]      5th order rosette aberration (�)
    aberrations.obj_lens_phi_54 = 0.00;                        % [phi_R5]	Azimuthal angle of 5th order rosette aberration (deg)
    aberrations.obj_lens_c_56 = 0.00;                          % [A5]      6-fold astigmatism (�)
    aberrations.obj_lens_phi_56 = 0.00;                        % [phi_A5]	Azimuthal angle of 6-fold astigmatism (deg)
end