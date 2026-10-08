function run_performance_50million_20261004
% Same full-array transport as the accepted benchmark; isolated extra scale.
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'performance_50million_20261004');
runtime=fullfile(out,'final_runtime');
assert(exist(runtime,'dir')==7,'Prepare the isolated runtime before starting.');
assert(~exist(fullfile(out,'run_started.txt'),'file'),'This run directory has already been used.');
cd(out); addpath(runtime,'-begin');
diary(fullfile(out,'run_diary.log')); cleanup=onCleanup(@() diary('off'));
write_text(fullfile(out,'matlab_pid.txt'),sprintf('%d',feature('getpid')));
write_text(fullfile(out,'run_started.txt'),'Dedicated 50-million full-array benchmark started.');
copyfile([mfilename('fullpath') '.m'],fullfile(out,'benchmark_driver.m'));
paths=cell(0,3); files=dir(fullfile(runtime,'*.m'));
for k=1:numel(files)
    name=files(k).name(1:end-2); preferred=which(name);
    assert(strcmpi(preferred,fullfile(runtime,files(k).name)),'Unexpected runtime path.');
    paths(end+1,:)={name,preferred,strjoin(which(name,'-all'),' | ')}; %#ok<AGROW>
end
writetable(cell2table(paths,'VariableNames',{'function_name','preferred_path','all_paths'}),fullfile(out,'matlab_function_paths.csv'));
config=struct('N',50000000,'warmups',2,'formal_repeats',5,'seed',20261004, ...
    'I_beam_A',0,'ee_method','self_consistent','ee_steps',4000,'phase','midpoint', ...
    'physics','nonrelativistic','path_mode','original_full_array', ...
    'memory_guard_factor',1.15,'physical_reserve_bytes',16*2^30, ...
    'watchdog','available physical RAM <8GiB for 5 consecutive 2s samples stops this MATLAB PID only', ...
    'output_mode','MAT v7.3 pidx/x/y per repeat; final retained');
write_text(fullfile(out,'benchmark_config.json'),jsonencode(config,'PrettyPrint',true));
[~,memsys]=memory;
hw=struct('MATLAB',version,'release',version('-release'),'computer',computer, ...
    'maxNumCompThreads',maxNumCompThreads,'OMP_NUM_THREADS',getenv('OMP_NUM_THREADS'), ...
    'MKL_NUM_THREADS',getenv('MKL_NUM_THREADS'),'available_RAM_bytes',memsys.PhysicalMemory.Available, ...
    'total_RAM_bytes',memsys.PhysicalMemory.Total,'client_date','2026-10-04','timezone','Asia/Shanghai');
write_text(fullfile(out,'matlab_hardware_readback.json'),jsonencode(hw,'PrettyPrint',true));
try
    validate_runtime(out);
    old=readtable(fullfile(root,'performance_final_20261004','tmm_timing_raw.csv'));
    large=old(old.N==10000000,:);
    increment_per_N=max((large.max_sampled_memory_bytes-large.baseline_memory_bytes)./large.N);
    records=struct([]);
    for rep=-config.warmups+1:config.formal_repeats
        [base,available]=memory_bytes;
        estimated_increment=config.memory_guard_factor*increment_per_N*config.N;
        gate=struct('repeat',rep,'baseline_memory_bytes',base,'available_RAM_bytes',available, ...
            'estimated_increment_bytes',estimated_increment,'reserve_bytes',config.physical_reserve_bytes, ...
            'pass',estimated_increment<=available-config.physical_reserve_bytes);
        write_text(fullfile(out,'latest_memory_gate.json'),jsonencode(gate,'PrettyPrint',true));
        assert(gate.pass,'Available RAM cannot accommodate the measured-growth estimate plus 16GiB reserve.');
        update_status(out,'running',rep,numel(records),'initializing');
        fprintf('FIFTY_MILLION_START rep=%d estimate_increment=%.2fGiB available=%.2fGiB\n',rep,estimated_increment/2^30,available/2^30);
        row=one_repeat(out,config.N,rep,config.seed);
        records=[records;row]; %#ok<AGROW>
        writetable(struct2table(records,'AsArray',true),fullfile(out,'tmm_timing_raw.csv'));
        update_status(out,'running',rep,numel(records),'repeat_complete');
        fprintf('FIFTY_MILLION_REPEAT rep=%d core=%.6fs full=%.6fs sampled=%.2fGiB\n',rep,row.core_s,row.full_wall_s,row.max_sampled_memory_bytes/2^30);
    end
    summarize(out,root);
    update_status(out,'complete',5,numel(records),'complete');
    write_text(fullfile(out,'PERFORMANCE_50M_COMPLETE.txt'),'Two warmups and five formal repeats completed, all outputs finite.');
    fprintf('FIFTY_MILLION_ALL_COMPLETE\n');
catch ME
    write_text(fullfile(out,'matlab_failure.txt'),getReport(ME,'extended','hyperlinks','off'));
    update_status(out,'failed',NaN,0,ME.message);
    rethrow(ME);
end
clear cleanup;
end

function validate_runtime(out)
[V,p]=load_comsol_initial_csv(fullfile(out,'inputs','tmm_canonical_initial.csv'));
assert(isequal(sort(p.id(:)),(1:50000)'));
p.epsilon_ee=zeros(1,50000); [sys,err]=case06_config;
[xy,~,ee,~]=transport(p,V,sys,err,50000);
m=metrics(xy);
assert(~ee.hit_radius_floor && ee.scale==1 && all(ee.delta_R==0));
ref=sortrows(readtable(fullfile(out,'reference_results','tmm_image_with_pidx.csv')),'pidx');
assert(isequal(ref.pidx,p.id(:)));
delta=max(abs(xy-ref{:,2:3}),[],'all')*1e9;
assert(delta<1e-8 && abs(m.fw99-39.0551665557053)<1e-8);
write_text(fullfile(out,'validation_summary.json'),jsonencode(struct('pass',true, ...
    'Case06_coordinate_max_component_difference_nm',delta,'Case06_FW99_nm',m.fw99, ...
    'runtime_scope','approved merged EB0/E1/ee; old B; midpoint phase; I0'), 'PrettyPrint',true));
fprintf('FIFTY_MILLION_VALIDATION_PASS\n');
end

function [sys,err]=case06_config
[sys,err]=table_iii_tmm_config('06_b_field_gradient_edge','legacy');
sys.func_B_vec=@(X,Y,Z) [-3.*(X-1e-6);-3.*Y;1.2+6.*(Z-0.5e-3)];
sys.physics_model='nonrelativistic';
if isfield(sys,'b_gradient'), sys=rmfield(sys,'b_gradient'); end
if isfield(sys,'nonrel_phase'), sys=rmfield(sys,'nonrel_phase'); end
err.I_beam=0; err.ee_method='self_consistent'; err.ee_steps=4000; err.ee_quiet=true;
assert(~isfield(err,'ee_reference_profile'));
end

function row=one_repeat(out,N,rep,seed)
[baseline,~]=memory_bytes; total=tic; load_s=0;
cfg=struct('load_from_file',false,'pattern_mode','spot','spot_radius',18e-9, ...
    'spot_center_x',0,'spot_center_y',0,'U0_max',0.5,'alpha_max_deg',10.7);
rng(seed,'twister'); a=tic; [V,p]=init_particles(N,cfg); p.id=1:N; init_s=toc(a);
p.epsilon_ee=zeros(1,N); [sys,err]=case06_config;
update_status(out,'running',rep,rep+1,'transport');
[xy,t,ee,sampled]=transport(p,V,sys,err,N);
assert(ee.scale==1 && ~ee.hit_radius_floor && all(ee.delta_R==0));
update_status(out,'running',rep,rep+1,'statistics');
a=tic; m=metrics(xy); statistics_s=toc(a);
update_status(out,'running',rep,rep+1,'output');
a=tic; pidx=p.id(:); x_m=xy(:,1); y_m=xy(:,2); %#ok<NASGU>
save(fullfile(out,sprintf('N%d_latest_image.mat',N)),'pidx','x_m','y_m','-v7.3'); output_s=toc(a);
full_wall=toc(total);
row=struct('N',N,'repeat',rep,'is_warmup',rep<=0,'condition','Case06 axis+1um NR I0', ...
    'input_mode','original generator; new unpaired scaling sample','seed',seed,'I_beam_A',0, ...
    'ee_steps',4000,'path_mode','original_full_array','load_s',load_s,'init_s',init_s, ...
    'common_physics_s',t.common,'operator_build_s',t.operators,'matrix_mapping_s',t.mapping, ...
    'core_s',t.common+t.operators+t.mapping,'core_cpu_s',t.cpu,'statistics_s',statistics_s, ...
    'output_s',output_s,'full_wall_s',full_wall,'cx_nm',m.cx,'cy_nm',m.cy, ...
    'RMS_radius_nm',m.rms,'FW99_nm',m.fw99,'baseline_memory_bytes',baseline, ...
    'max_sampled_memory_bytes',max(sampled),'phase_method','midpoint', ...
    'ee_scale',ee.scale,'ee_floor_hit',ee.hit_radius_floor);
if rep==5
    sample_index=unique(round(linspace(1,N,1000)));
    V_sample=V(:,:,sample_index); p_sample=struct;
    fields=fieldnames(p);
    for k=1:numel(fields), p_sample.(fields{k})=p.(fields{k})(sample_index); end
    xy_sample=xy(sample_index,:); %#ok<NASGU>
    save(fullfile(out,'validation_sample.mat'),'sample_index','V_sample','p_sample','xy_sample');
    [check,~,~,~]=transport(p_sample,V_sample,sys,err,numel(sample_index));
    validation_difference=max(abs(check-xy_sample),[],'all')*1e9;
    assert(validation_difference<1e-8,'Selected-particle full-array mapping check failed.');
    write_text(fullfile(out,'sample_mapping_validation.json'),jsonencode(struct( ...
        'samples',numel(sample_index),'maximum_component_difference_nm',validation_difference, ...
        'scope','I0 native ee scale exactly1; selected actual input replay; outside reported full_wall'), 'PrettyPrint',true));
end
end

function [xy,t,ee,sampled]=transport(p,V,sys,err,N)
cpu=cputime; sampled=zeros(1,4); sampled(1)=memory_bytes;
a=tic; phys=calc_common_physics(p,V,sys); t.common=toc(a); sampled(2)=memory_bytes;
assert(~phys.is_relativistic && strcmp(phys.nonrel.phase_diagnostics.method,'midpoint'));
a=tic;
M_spread=calc_M_spread(p,phys,N);
M_Br1=calc_M_Br1(p,phys,sys,err,N);
M_B1=calc_M_B1(p,phys,sys,err,N);
M_EBrho=calc_M_EB_angle(p,phys,sys,err,N);
M_EB0=calc_M_EB0(p,phys,sys,err,N);
[M_ee,ee]=calc_M_ee(p,phys,sys,err,N);
M_E1=calc_M_E1(p,phys,sys,err,N);
t.operators=toc(a); sampled(3)=memory_bytes;
a=tic;
M_temp1=pagemtimes(M_EB0,M_spread);
M_temp2=pagemtimes(M_Br1,M_temp1);
M_temp3=pagemtimes(M_E1,M_temp2);
M_temp4=pagemtimes(M_B1,M_temp3);
M_temp5=pagemtimes(M_EBrho,M_temp4);
M_total=pagemtimes(M_ee,M_temp5);
V_final=pagemtimes(M_total,V); Pos_final=squeeze(V_final);
xy=Pos_final(1:2,:).'; t.mapping=toc(a); sampled(4)=memory_bytes;
t.cpu=cputime-cpu;
end

function summarize(out,root)
T=readtable(fullfile(out,'tmm_timing_raw.csv')); x=T(T.repeat>0,:);
assert(height(x)==5 && sum(T.repeat<=0)==2 && all(T.N==50000000));
assert(max(T.FW99_nm)-min(T.FW99_nm)<1e-10,'Repeat outputs differ with the same seed.');
row=struct('N',50000000,'warmups',2,'formal_repeats',5, ...
    'core_median_s',median(x.core_s),'core_min_s',min(x.core_s),'core_max_s',max(x.core_s), ...
    'core_IQR_s',qlinear(x.core_s,0.75)-qlinear(x.core_s,0.25), ...
    'load_median_s',median(x.load_s),'init_median_s',median(x.init_s), ...
    'common_median_s',median(x.common_physics_s),'operators_median_s',median(x.operator_build_s), ...
    'mapping_median_s',median(x.matrix_mapping_s),'statistics_median_s',median(x.statistics_s), ...
    'output_median_s',median(x.output_s),'full_wall_median_s',median(x.full_wall_s), ...
    'max_stage_sampled_memory_bytes',max(x.max_sampled_memory_bytes));
S=struct2table(row,'AsArray',true); writetable(S,fullfile(out,'tmm_timing_summary.csv'));
previous=readtable(fullfile(root,'performance_final_20261004','tmm_timing_summary.csv'));
combined=[previous;S]; writetable(combined,fullfile(out,'tmm_timing_summary_with_50million.csv'));
all_raw=[readtable(fullfile(root,'performance_final_20261004','tmm_timing_raw.csv'));T];
writetable(all_raw,fullfile(out,'tmm_timing_raw_with_50million.csv'));
f=figure('Visible','off','Color','w','Position',[100,100,820,520]);
loglog(combined.N,combined.core_median_s,'o-','LineWidth',1.4); hold on;
errorbar(combined.N,combined.core_median_s,combined.core_median_s-combined.core_min_s,combined.core_max_s-combined.core_median_s,'LineStyle','none');
xlabel('Number of particles'); ylabel('NR-TMM transport core (s)'); grid on;
print(f,fullfile(out,'tmm_scaling_with_50million.png'),'-dpng','-r220'); close(f);
fid=fopen(fullfile(out,'performance_report.md'),'w','n','UTF-8');
fprintf(fid,'# Measured 50-Million Full-Array NR-TMM\n\nTwo warmups and five formal repeats, unchanged frozen runtime and full cascade.\n');
fprintf(fid,'Core median %.9g s; min %.9g s; max %.9g s; IQR %.9g s.\n',row.core_median_s,row.core_min_s,row.core_max_s,row.core_IQR_s);
fprintf(fid,'Full wall median %.9g s; stage-boundary memory maximum %.6g GiB (not continuous peak).\n',row.full_wall_median_s,row.max_stage_sampled_memory_bytes/2^30);
fprintf(fid,'FW99 %.15g nm, cx %.15g nm, cy %.15g nm; all 50 million coordinates finite.\n',T.FW99_nm(end),T.cx_nm(end),T.cy_nm(end));
fprintf(fid,'Seed20261004 original spot generator, new unpaired scale sample; not a COMSOL shared input. Native ee4000/I0, no fast return, GPU, batching or thread override.\n');
fprintf(fid,'Guard reduced from the earlier conservative factor1.6/85%% criterion to factor1.15 measured-growth increment plus16GiB reserve; external2s monitor samples RAM/process memory.\n');
fprintf(fid,'Full wall includes init, transport, statistics, output and progress metadata; excludes MATLAB startup, preflight acceptance and final selected-particle validation.\n');
fprintf(fid,'Previous 50-million point was not attempted; this new actual success supersedes only that scale status. No new COMSOL timing or memory-wall evidence.\n'); fclose(fid);
end

function m=metrics(xy)
assert(all(isfinite(xy),'all'),'Nonfinite output.'); c=mean(xy,1);
r=hypot(xy(:,1)-c(1),xy(:,2)-c(2))*1e9;
m=struct('cx',c(1)*1e9,'cy',c(2)*1e9,'rms',sqrt(mean(r.^2)),'fw99',2*qlinear(r,0.99));
end
function q=qlinear(x,p)
x=sort(x(:)); at=1+(numel(x)-1)*p; lo=floor(at); hi=ceil(at);
q=x(lo)+(at-lo)*(x(hi)-x(lo));
end
function [bytes,available]=memory_bytes
[user,sys]=memory; bytes=user.MemUsedMATLAB; available=sys.PhysicalMemory.Available;
end
function update_status(out,status,rep,completed,stage)
write_text(fullfile(out,'execution_status.json'),jsonencode(struct('status',status, ...
    'N',50000000,'repeat',rep,'completed_repeats',completed,'stage',stage),'PrettyPrint',true));
end
function write_text(file,txt)
fid=fopen(file,'w','n','UTF-8'); assert(fid>=0); fprintf(fid,'%s\n',txt); fclose(fid);
end
