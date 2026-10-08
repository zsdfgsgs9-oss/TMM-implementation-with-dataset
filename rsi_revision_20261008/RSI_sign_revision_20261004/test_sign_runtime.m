function test_sign_runtime()
here=fileparts(mfilename('fullpath')); root=fileparts(here);
addpath(fullfile(here,'runtime'),'-begin'); clear calc_M_Br1
[sys,err]=table_iii_tmm_config('06_b_field_gradient_axis_1um');
assert(~isfield(sys,'b_gradient') && ~isfield(sys,'nonrel_phase'));
E=sys.func_E_vec(0,0); assert(E(3)==-sys.U/sys.L);
p=struct('U0',.5,'alpha',.1,'phi',0,'epsilon_ee',0);
v=reshape([0;0;1],[3 1 1]); phys=calc_common_physics(p,v,sys);
assert(phys.nonrel.a_z>0);
M=calc_M_Br1(p,phys,sys,err,1);
assert(M(1,2)>0 && M(2,1)<0 && M(1,2)==-M(2,1));
B=[0;0;1.2]; assert(isequal(cross([2;3;-1],B)/1.2^2,[2.5;-5/3;0]));
assert((cross([0;0;-sys.U/sys.L],[-3e-6;0;1.2]))'*[0;1;0]>0);
[sys0,err0]=table_iii_tmm_config('01_intrinsic_error');
out=run_transfer_matrix_model(p,v,sys0,err0,1);
assert(~out.diagnostics.b_gradient.enabled);
% Match the existing display selection without altering the launched ensemble.
stream=RandStream('mt19937ar','Seed',2026100408);
pick=sort(randperm(stream,50000,3000));
writetable(table(pick(:),'VariableNames',{'pidx'}),fullfile(here,'stripe_display_ids.csv'));
initial=fullfile(root,'RSI_revision_simulation_package_20260901', ...
    '07_source_inputs','intrinsic_pair','01_intrinsic_error_initial.csv');
[v,p]=load_comsol_initial_csv(initial); p.epsilon_ee=zeros(size(p.U0));
[sys,err]=table_iii_tmm_config('06_b_field_gradient_axis_1um');
before_br=fullfile(here,'old_br'); seconds=zeros(10,2);
for rep=-2:10
    for j=1:2
        cfg=sys;
        if j==1
            E=sys.U/sys.L; cfg.func_E_vec=@(X,Y)[zeros(size(X));zeros(size(Y));E*ones(size(X))];
            addpath(before_br,'-begin');
        else
            if contains(path,before_br), rmpath(before_br); end
        end
        clear calc_M_Br1
        start=tic; evalc('out=run_transfer_matrix_model(p,v,cfg,err,50000);');
        elapsed=toc(start);
        if rep>0, seconds(rep,j)=elapsed; end
    end
end
writetable(table((1:10)',seconds(:,1),seconds(:,2), ...
    'VariableNames',{'repeat','before_call_s','signed_call_s'}),fullfile(here,'local_timing_control.csv'));
fprintf('Signed operator and direction tests passed; local before/after medians %.6g / %.6g s.\n',median(seconds));
end
