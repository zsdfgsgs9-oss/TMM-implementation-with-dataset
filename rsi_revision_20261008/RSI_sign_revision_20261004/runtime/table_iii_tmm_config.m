function [sys, err] = table_iii_tmm_config(case_name, variant)
% Published near-axis operators; E_vec is the physical Cartesian field.
if nargin > 1 && ~any(string(variant)==["published","legacy"])
    error('This runtime contains only the published, separate B operators.');
end
sys=struct('B0',1.2,'L',1e-3,'U',12.67e3,'N_turn',1, ...
    'physics_model','nonrelativistic');
E=sys.U/sys.L;
sys.func_E_vec=@(X,Y)[zeros(size(X));zeros(size(Y));-E*ones(size(X))];
sys.func_B_vec=@(X,Y,Z)[zeros(size(X));zeros(size(Y));sys.B0*ones(size(Z))];
err=struct('G',0,'delta_E',0,'delta_B',0,'z1',0,'z2',sys.L, ...
    'I_beam',0,'R0',18e-9,'ee_method','self_consistent','ee_steps',4000);
err.func_dEz=@(x,y,z)zeros(size(x)); % Positive accelerating-field perturbation.
switch string(case_name)
    case "01_intrinsic_error"
    case "02_em_field_angle_error"
        sys.func_E_vec=@(X,Y)[.8e-3*E*ones(size(X)); ...
            .6e-3*E*ones(size(Y));-E*ones(size(X))];
    case "03_e_field_zero_order"
        err.delta_E=1.267e4;
    case "04_b_field_zero_order"
        err.delta_B=1.2e-3;
    case {"05_b_field_gradient","06_b_field_gradient_axis_1um","06_b_field_gradient_edge"}
        err.G=6; axis_x=0;
        if startsWith(string(case_name),"06"), axis_x=1e-6; end
        sys.func_B_vec=@(X,Y,Z)[-3*(X-axis_x);-3*Y;sys.B0+6*(Z-.5*sys.L)];
    case "07_e_field_gradient"
        err.func_dEz=@(x,y,z)6.335e7*z.*ones(size(x));
    otherwise
        error('Unknown published case: %s',case_name);
end
end
