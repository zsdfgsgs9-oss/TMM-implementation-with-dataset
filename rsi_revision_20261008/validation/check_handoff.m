function check_handoff(root, output_dir)
% Recheck retained data without changing scientific result files or solving FEM.
if nargin < 1, root = fileparts(fileparts(mfilename('fullpath'))); end
if nargin < 2, output_dir = tempname; end
assert(~strcmp(root,output_dir),'Use a separate output directory.');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
signed=fullfile(root,'RSI_sign_revision_20261004');
runtime=fullfile(signed,'runtime');
previous_path=path; cleanup=onCleanup(@()path(previous_path)); %#ok<NASGU>
addpath(runtime,'-begin');
clear run_transfer_matrix_model calc_common_physics calc_M_Br1 calc_M_E1 calc_M_EB0 table_iii_tmm_config
assert(strcmp(which('run_transfer_matrix_model'),fullfile(runtime,'run_transfer_matrix_model.m')));
[sys,err]=table_iii_tmm_config('06_b_field_gradient_axis_1um');
assert(strcmp(sys.physics_model,'nonrelativistic'));
assert(~isfield(sys,'b_gradient') && ~isfield(sys,'nonrel_phase'));
E=sys.func_E_vec(0,0); assert(E(3)==-sys.U/sys.L);
p1=struct('U0',.5,'alpha',.1,'phi',0,'epsilon_ee',0);
v1=reshape([0;0;1],[3 1 1]); phys=calc_common_physics(p1,v1,sys);
assert(phys.nonrel.a_z>0);
M=calc_M_Br1(p1,phys,sys,err,1);
assert(M(1,2)>0 && M(2,1)<0 && M(1,2)==-M(2,1));
assert(isequal(cross([2;3;-1],[0;0;1.2])/1.2^2,[2.5;-5/3;0]));
result=struct('operator_direction_tests_pass',true,'MATLAB_version',version, ...
    'COMSOL_solved',false,'performance_benchmarked',false);
initial=fullfile(root,'RSI_revision_simulation_package_20260901', ...
    '07_source_inputs','intrinsic_pair','01_intrinsic_error_initial.csv');
[v,p]=load_comsol_initial_csv(initial); p.epsilon_ee=zeros(size(p.U0));
names=["01_intrinsic_error","02_em_field_angle_error","03_e_field_zero_order", ...
    "04_b_field_zero_order","05_b_field_gradient","06_b_field_gradient_axis_1um", ...
    "07_e_field_gradient"];
replay_error=zeros(1,7);
for k=1:7
    pk=p; vk=v;
    if k==6
        input6=fullfile(root,'RSI非空间电荷执行结果_20261001','02_旧矩阵回归', ...
            'Case06_axis_1um_20261002','initial_raw.csv');
        [vk,pk]=load_comsol_initial_csv(input6); pk.epsilon_ee=zeros(size(pk.U0));
    end
    [sys,err]=table_iii_tmm_config(names(k));
    out=run_transfer_matrix_model(pk,vk,sys,err,50000);
    assert(~out.diagnostics.b_gradient.enabled);
    [ids,ix]=sort(pk.id(:)); assert(isequal(ids,(1:50000)'));
    retained=sortrows(readtable(fullfile(signed,sprintf('case%02d_coordinates.csv',k))),'pidx');
    assert(isequal(retained.pidx,ids));
    xy=[out.X_final(ix).' out.Y_final(ix).']*1e9;
    replay_error(k)=max(abs(xy-[retained.tmm_x_nm retained.tmm_y_nm]),[],'all');
    assert(replay_error(k)<1e-7,'Signed case %d replay differs.',k);
    fprintf('Signed case %d replay max error %.9g nm\n',k,replay_error(k));
end
result.seven_case_replay_max_error_nm=replay_error;

source=fullfile(root,'Fig10最小工况重跑结果_20261002');
retained=readtable(fullfile(signed,'stripe_paired_coordinates.csv'),'TextType','string');
stripe_names={'11_center','13_axis_minus_1um'};
input_names={'11_system_stripe_center','13_system_stripe_1um'};
offsets=[0 1e-6]; stripe_error=zeros(1,2);
sys=struct('B0',1.2,'L',1e-3,'U',12.67e3,'N_turn',1,'physics_model','nonrelativistic');
E=sys.U/sys.L; G=6;
sys.func_E_vec=@(X,Y)[.8e-3*E*ones(size(X));.6e-3*E*ones(size(Y));-E*ones(size(X))];
sys.func_B_vec=@(X,Y,Z)[-.5*G*X;-.5*G*Y;sys.B0+G*(Z-.5*sys.L)];
err=struct('G',G,'delta_E',1.267e4,'delta_B',1.2e-3,'z1',0, ...
    'z2',sys.L,'I_beam',0,'R0',18e-9,'ee_method','gauss');
err.func_dEz=@(x,y,z)6.335e7*z.*ones(size(x));
for k=1:2
    raw=readmatrix(fullfile(source,'inputs',[input_names{k} '_initial.csv']),'FileType','text');
    raw=raw(all(isfinite(raw(:,1:7)),2),1:7);
    assert(size(raw,1)==50000 && isequal(round(raw(:,1)),(1:50000)'));
    pk=struct('id',round(raw(:,1)).','U0',raw(:,4).'/1.602e-19, ...
        'alpha',atan2(hypot(raw(:,5),raw(:,6)),raw(:,7)).', ...
        'phi',atan2(raw(:,6),raw(:,5)).','epsilon_ee',zeros(1,50000));
    vk=zeros(3,1,50000);
    vk(1,1,:)=raw(:,2).'*1e-6+offsets(k);
    vk(2,1,:)=raw(:,3).'*1e-6; vk(3,1,:)=1;
    out=run_transfer_matrix_model(pk,vk,sys,err,50000);
    assert(~out.diagnostics.b_gradient.enabled);
    xy=[out.X_final(:)-offsets(k),out.Y_final(:)]*1e9;
    pair=sortrows(retained(retained.case_name==string(stripe_names{k}),:),'pidx');
    stripe_error(k)=max(abs(xy-[pair.tmm_x_nm pair.tmm_y_nm]),[],'all');
    assert(stripe_error(k)<1e-7,'Stripe replay differs.');
end
result.stripe_replay_max_error_nm=stripe_error;

quad=fullfile(root,'RSI最终版报告与汇总数据_精简版_20261004','02_双场口径rtol1e7');
origins={'entrance','midpoint'}; cases={'H0_control','H0p5','H1','H2'};
xy=zeros(50000,2,4,2);
for g=1:2
    for k=1:4
        a=sortrows(readtable(fullfile(quad,origins{g},cases{k},'image_particles.csv')),'pidx');
        assert(isequal(a.pidx,(1:50000)'));
        xy(:,:,k,g)=[a.x_m a.y_m]*1e9;
    end
end
% Use the original Twister seed and the same shared index draw for all cases.
rng(20261004,'twister'); nboot=1000; bootstrap=zeros(nboot,4,2);
for b=1:nboot
    ix=randi(50000,50000,1);
    for g=1:2
        f0=radial_fw99(xy(ix,:,1,g));
        for k=2:4
            bootstrap(b,k,g)=radial_fw99(xy(ix,:,k,g))-f0;
        end
    end
    if mod(b,100)==0, fprintf('Quadratic paired bootstrap %d/%d\n',b,nboot); end
end
ci=readtable(fullfile(quad,'paired_bootstrap_CI95.csv'),'TextType','string');
intervals=struct([]); max_ci_error=0;
for g=1:2
    for k=1:4
        row=ci(ci.field_reference==string(origins{g}) & ci.case_name==string(cases{k}),:);
        assert(height(row)==1 && row.seed==20261004 && row.bootstrap_repeats==nboot);
        limits=qlinear(bootstrap(:,k,g),[.025 .975]);
        difference=max(abs(limits-[row.H_H0_CI95_low_nm row.H_H0_CI95_high_nm]));
        max_ci_error=max(max_ci_error,difference);
        assert(difference<1e-9,'Quadratic paired CI differs from retained interval.');
        intervals=[intervals;struct('field_reference',origins{g},'case_name',cases{k}, ...
            'ci95_low_nm',limits(1),'ci95_high_nm',limits(2),'max_abs_error_nm',difference)]; %#ok<AGROW>
    end
end
result.quadratic_bootstrap_repeats=nboot;
result.quadratic_bootstrap_seed=20261004;
result.quadratic_H_minus_H0_intervals=intervals;
result.quadratic_CI_max_abs_error_nm=max_ci_error;
result.other_5000_repeat_intervals_rerun=false;
fid=fopen(fullfile(output_dir,'matlab_validation.json'),'w','n','UTF-8');
assert(fid>=0); closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(result,'PrettyPrint',true));
fprintf('Signed replay, direction tests and quadratic CI verification passed.\n');
end

function value=radial_fw99(xy)
c=mean(xy,1); r=hypot(xy(:,1)-c(1),xy(:,2)-c(2));
value=2*qlinear(r,.99);
end

function q=qlinear(x,p)
x=sort(x(:)); pos=1+(numel(x)-1)*p;
q=reshape(x(floor(pos)),size(p))+(pos-floor(pos)).*reshape(x(ceil(pos))-x(floor(pos)),size(p));
end
