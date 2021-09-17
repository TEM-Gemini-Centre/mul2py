function [input_multem] = set_aberrations(input_multem, aberrations)
    %% CL aberrations
    input_multem.cond_lens_c_10 = aberrations.cond_lens_c_10;
    input_multem.cond_lens_c_12 = aberrations.cond_lens_c_12;
    input_multem.cond_lens_phi_12 = aberrations.cond_lens_phi_12;

    input_multem.cond_lens_c_21 = aberrations.cond_lens_c_21;
    input_multem.cond_lens_phi_21 = aberrations.cond_lens_phi_21;
    input_multem.cond_lens_c_23 = aberrations.cond_lens_c_23;
    input_multem.cond_lens_phi_23 = aberrations.cond_lens_phi_23;

    input_multem.cond_lens_c_30 = aberrations.cond_lens_c_30;
    input_multem.cond_lens_c_32 = aberrations.cond_lens_c_32;
    input_multem.cond_lens_phi_32 = aberrations.cond_lens_phi_32;
    input_multem.cond_lens_c_34 = aberrations.cond_lens_c_34;
    input_multem.cond_lens_phi_34 = aberrations.cond_lens_phi_34;

    input_multem.cond_lens_c_41 = aberrations.cond_lens_c_41;
    input_multem.cond_lens_phi_41 = aberrations.cond_lens_phi_41;
    input_multem.cond_lens_c_43 = aberrations.cond_lens_c_43;
    input_multem.cond_lens_phi_43 = aberrations.cond_lens_phi_43;
    input_multem.cond_lens_c_45 = aberrations.cond_lens_c_45;
    input_multem.cond_lens_phi_45 = aberrations.cond_lens_phi_45;

    input_multem.cond_lens_c_50 = aberrations.cond_lens_c_50;
    input_multem.cond_lens_c_52 = aberrations.cond_lens_c_52;
    input_multem.cond_lens_phi_52 = aberrations.cond_lens_phi_52;
    input_multem.cond_lens_c_54 = aberrations.cond_lens_c_54;
    input_multem.cond_lens_phi_54 = aberrations.cond_lens_phi_54;
    input_multem.cond_lens_c_56 = aberrations.cond_lens_c_56;
    input_multem.cond_lens_phi_56 = aberrations.cond_lens_phi_56;

    %% OL aberrations
    input_multem.obj_lens_c_10 = aberrations.obj_lens_c_10;
    input_multem.obj_lens_c_12 = aberrations.obj_lens_c_12;
    input_multem.obj_lens_phi_12 = aberrations.obj_lens_phi_12;

    input_multem.obj_lens_c_21 = aberrations.obj_lens_c_21;
    input_multem.obj_lens_phi_21 = aberrations.obj_lens_phi_21;
    input_multem.obj_lens_c_23 = aberrations.obj_lens_c_23;
    input_multem.obj_lens_phi_23 = aberrations.obj_lens_phi_23;

    input_multem.obj_lens_c_30 = aberrations.obj_lens_c_30;
    input_multem.obj_lens_c_32 = aberrations.obj_lens_c_32;
    input_multem.obj_lens_phi_32 = aberrations.obj_lens_phi_32;
    input_multem.obj_lens_c_34 = aberrations.obj_lens_c_34;
    input_multem.obj_lens_phi_34 = aberrations.obj_lens_phi_34;

    input_multem.obj_lens_c_41 = aberrations.obj_lens_c_41;
    input_multem.obj_lens_phi_41 = aberrations.obj_lens_phi_41;
    input_multem.obj_lens_c_43 = aberrations.obj_lens_c_43;
    input_multem.obj_lens_phi_43 = aberrations.obj_lens_phi_43;
    input_multem.obj_lens_c_45 = aberrations.obj_lens_c_45;
    input_multem.obj_lens_phi_45 = aberrations.obj_lens_phi_45;

    input_multem.obj_lens_c_50 = aberrations.obj_lens_c_50;
    input_multem.obj_lens_c_52 = aberrations.obj_lens_c_52;
    input_multem.obj_lens_phi_52 = aberrations.obj_lens_phi_52;
    input_multem.obj_lens_c_54 = aberrations.obj_lens_c_54;
    input_multem.obj_lens_phi_54 = aberrations.obj_lens_phi_54;
    input_multem.obj_lens_c_56 = aberrations.obj_lens_c_56;
    input_multem.obj_lens_phi_56 = aberrations.obj_lens_phi_56;
    
end
