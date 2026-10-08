function [sys, err] = table_iii_tmm_config(case_name, variant)
%TABLE_III_TMM_CONFIG Reproducible Table III nonrelativistic TMM settings.

    if nargin < 2 || strlength(string(variant)) == 0
        variant = "corrected";
    end
    case_name = string(case_name);
    variant = lower(string(variant));

    sys.B0 = 1.2;
    sys.L = 1e-3;
    sys.U = 12.67e3;
    sys.N_turn = 1;
    sys.physics_model = "nonrelativistic";
    E = sys.U / sys.L;
    sys.func_E_vec = @(X, Y) [zeros(size(X)); zeros(size(Y)); ...
        E .* ones(size(X))];
    sys.func_B_vec = @(X, Y, Z) [zeros(size(X)); zeros(size(Y)); ...
        sys.B0 .* ones(size(Z))];
    err = zero_errors(sys.L);

    switch case_name
        case "01_intrinsic_error"
        case "02_em_field_angle_error"
            sys.func_E_vec = @(X, Y) [0.8e-3 .* E .* ones(size(X)); ...
                0.6e-3 .* E .* ones(size(Y)); E .* ones(size(X))];
        case "03_e_field_zero_order"
            err.delta_E = 1.267e4;
        case "04_b_field_zero_order"
            err.delta_B = 1.2e-3;
        case {"05_b_field_gradient", "06_b_field_gradient_edge"}
            err.G = 6;
            axis_center = [0, 0];
            if case_name == "06_b_field_gradient_edge"
                axis_center = [1e-3, 0];
            end
            sys.func_B_vec = @(X, Y, Z) [-3 .* (X - axis_center(1)); ...
                -3 .* (Y - axis_center(2)); ...
                sys.B0 + 6 .* (Z - 0.5 .* sys.L)];
            sys.nonrel_phase.method = "linear";
            sys.nonrel_phase.G_Tpm = 6;
            if variant ~= "legacy" && variant ~= "phase_only"
                sys.b_gradient.enabled = true;
                sys.b_gradient.G_Tpm = 6;
                sys.b_gradient.axis_center_m = axis_center;
                sys.b_gradient.integration_steps = 128;
                sys.b_gradient.physical_E_vec = [0; 0; -E];
                sys.b_gradient.include_axis_offset = true;
                sys.b_gradient.include_gradient_force = true;
                sys.b_gradient.include_phase_variation = true;
                sys.b_gradient.include_transverse_e = true;
                if variant == "unified_no_axis"
                    sys.b_gradient.include_axis_offset = false;
                elseif variant == "unified_no_gradient_force"
                    sys.b_gradient.include_gradient_force = false;
                elseif variant == "unified_midpoint_phase"
                    sys.b_gradient.include_phase_variation = false;
                elseif variant == "unified_256"
                    sys.b_gradient.integration_steps = 256;
                elseif variant ~= "corrected" && variant ~= "unified"
                    error('Unknown B-gradient variant "%s".', variant);
                end
            elseif variant == "legacy"
                sys = rmfield(sys, 'nonrel_phase');
            end
        case "07_e_field_gradient"
            err.func_dEz = @(x, y, z) 6.335e7 .* z .* ones(size(x));
        otherwise
            error('Unknown Table III case "%s".', case_name);
    end
end

function err = zero_errors(L)
    err.G = 0;
    err.delta_E = 0;
    err.delta_B = 0;
    err.z1 = 0;
    err.z2 = L;
    err.func_dEz = @(x, y, z) zeros(size(x));
    err.I_beam = 0;
    err.R0 = 18e-9;
    err.ee_method = 'gauss';
end
