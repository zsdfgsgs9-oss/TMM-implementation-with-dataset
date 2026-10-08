function run_reviewer7_bootstrap(nboot)
% Paired percentile intervals for the manuscript's retained operators/planes.
if nargin < 1, nboot = 5000; end
here = fileparts(mfilename('fullpath')); root = fileparts(here);
archive = fullfile(root, 'RSI大修仿真统一归档_20260911');
third = fullfile(root, 'RSI第三版任务1-5成果归档_20260915', ...
    '02_任务1-3_MATLAB算子修正与回归');
old = readtable(fullfile(root, 'RSI第四版结论数据报告归档_20260916', ...
    '04_分项数据', 'D1_统计区间', 'table_iii_bootstrap.csv'), 'TextType', 'string');
rows = table();
for c = ["01_intrinsic_error", "02_em_field_angle_error", ...
        "03_e_field_zero_order"]
    a = old(old.case_name == c, :);
    rows = [rows; import_old(a)]; %#ok<AGROW>
end
case07=readtable(fullfile(root,'RSI_reviewer8_20261004', ...
    'case07_corrected_intervals.csv'),'TextType','string');
assert(all(case07.n_bootstrap==nboot) && all(case07.seed==2026100407));
rows=[rows;case07];
new = readtable(fullfile(root, 'RSI_EB0_recheck_20261003', ...
    'eb0_recheck_coordinates_with_pidx.csv'), 'TextType', 'string');
c = "04_b_field_zero_order";
t = sortrows(new(new.case_name == c, :), 'pidx');
f = load_comsol(fullfile(archive, '03_COMSOL完整结果', 'revision_comsol', ...
    'table_iii_nr_full_20260910', char(c), 'rel_off', ...
    'poincare_design_plane_with_pidx.csv'));
assert(isequal(t.pidx, f.id));
rows = [rows; one_case(here, "table_iii", c, ...
    {[t.x_nm t.y_nm], f.xy}, ["NR-TMM", "NR-full"], nboot, 2026100401)];

% Replay the frozen legacy B operator, without the alternative unified map.
addpath(fullfile(archive, '05_可复现脚本', 'root_matlab'), '-end');
addpath(fullfile(archive, '05_可复现脚本', 'revision_simulation'), '-end');
addpath(fullfile(third, 'code_snapshot'), '-begin');
addpath(fullfile(root, 'RSI_EB0_recheck_20261003', 'override'), '-begin');
clear calc_M_EB0
[v, p] = load_comsol_initial_csv(fullfile(root, ...
    'RSI_revision_simulation_package_20260901', '07_source_inputs', ...
    'intrinsic_pair', '01_intrinsic_error_initial.csv'));
p.epsilon_ee = zeros(1, numel(p.U0));
c = "05_b_field_gradient";
[sys, err] = table_iii_tmm_config(c, 'legacy');
out = run_transfer_matrix_model(p, v, sys, err, numel(p.U0));
[ids, order] = sort(p.id(:));
xy = [out.X_final(order).' out.Y_final(order).'] * 1e9;
f = load_comsol(fullfile(archive, '03_COMSOL完整结果', 'revision_comsol', ...
    'table_iii_nr_full_20260910', char(c), 'rel_off', ...
    'poincare_design_plane_with_pidx.csv'));
assert(isequal(ids, f.id));
m = beam_metrics(xy); mf = beam_metrics(f.xy);
assert(abs(m(7)-39.0541169682701) < 1e-6, 'Wrong legacy B operator.');
assert(abs(mf(7)-37.4389258239993) < 1e-6, 'Wrong Case 05 reference.');
writetable(table(ids, xy(:,1), xy(:,2), 'VariableNames', ...
    {'pidx','x_nm','y_nm'}), fullfile(here, 'case05_legacy_tmm_coordinates.csv'));
rows = [rows; one_case(here, "table_iii", c, {xy, f.xy}, ...
    ["NR-TMM", "NR-full"], nboot, 2026100402)];
base = fullfile(root,'RSI非空间电荷执行结果_20261001', ...
    '02_旧矩阵回归','Case06_axis_1um_20261002');
t = sortrows(readtable(fullfile(base,'tmm_image_with_pidx.csv')), 'pidx');
f = sortrows(readtable(fullfile(base,'fem_image_with_pidx.csv')), 'pidx');
assert(isequal(t.pidx,f.pidx) && numel(unique(t.pidx))==50000);
xy = {[t.x_m t.y_m]*1e9,[f.x_m f.y_m]*1e9};
m = beam_metrics(xy{1}); mf = beam_metrics(xy{2});
assert(abs(m(7)-39.0551665557053)<1e-6);
assert(abs(mf(7)-37.438920085442)<1e-6);
assert(abs(hypot(m(1)-mf(1),m(2)-mf(2))-2.95036935761173)<1e-6);
rows = [rows; one_case(here,"table_iii","06_b_field_gradient_axis_1um", ...
    xy,["NR-TMM","NR-full"],nboot,2026100403)];
rows = sortrows(rows,{'case_name','metric','model'});
writetable(rows, fullfile(here, 'table_iii_intervals.csv'));

hv = fullfile(root, 'HV0_12p67_to_150kV_results_20261003');
hvrows = table();
folders = ["U12p67kV", "U30kV", "U50kV", "U100kV", "U150kV"];
for k = 1:numel(folders)
    base = fullfile(hv, folders(k));
    t = sortrows(readtable(fullfile(base, 'TMM_at_TMM', 'tmm_image.csv')), 'pidx');
    nr = sortrows(readtable(fullfile(base, 'NR_at_TMM', 'image_particles.csv')), 'pidx');
    rel = sortrows(readtable(fullfile(base, 'REL_at_TMM', 'image_particles.csv')), 'pidx');
    ref = sortrows(readtable(fullfile(base, 'REL_own_focus', 'image_best_focus.csv')), 'pidx');
    assert(isequal(t.pidx, nr.pidx, rel.pidx, ref.pidx));
    xy = {[t.x_nm t.y_nm], [nr.x_m nr.y_m]*1e9, ...
        [rel.x_m rel.y_m]*1e9, [ref.x_m ref.y_m]*1e9};
    one = one_case(here, "voltage", folders(k), xy, ...
        ["NR-TMM", "NR-full", "REL-full", "REL-reference"], ...
        nboot, 2026100410+k);
    hvrows = [hvrows; one]; %#ok<AGROW>
    writetable(hvrows, fullfile(here, 'voltage_intervals.csv'));
end
fprintf('Completed retained-operator and fixed-plane bootstrap.\n');
end

function rows = one_case(here, dataset, c, xy, names, nboot, seed)
file = fullfile(here, sprintf('%s_%s_bootstrap.csv', dataset, c));
if isfile(file)
    rows = readtable(file, 'TextType', 'string');
    if all(rows.n_bootstrap == nboot) && all(rows.seed == seed), return; end
end
n = size(xy{1}, 1); assert(n == 50000);
nm = numel(xy); values = zeros(nm,7); boot = zeros(nboot,7,nm);
dc = zeros(nboot,1);
for j = 1:nm, values(j,:) = beam_metrics(xy{j}); end
fprintf('%s: %d models, %d paired replicates\n', c, nm, nboot);
for b = 1:nboot
    stream = RandStream('Threefry', 'Seed', seed); stream.Substream = b;
    idx = randi(stream,n,n,1);
    for j = 1:nm, boot(b,:,j) = beam_metrics(xy{j}(idx,:)); end
    dc(b) = hypot(boot(b,1,1)-boot(b,1,2), boot(b,2,1)-boot(b,2,2));
    if mod(b,1000)==0, fprintf('  %s: %d/%d\n', c,b,nboot); end
end
metrics = ["cx","cy","centroid_magnitude","RMS_radius","FW50","FW90","FW99"];
rows = table();
for j = 1:nm
    for q = 1:7
        rows = [rows; make_row(dataset,c,names(j),metrics(q),values(j,q), ...
            boot(:,q,j),nboot,seed)]; %#ok<AGROW>
    end
end
for pair = [1 2; 2 3; 1 3]'
    a = pair(1); bb = pair(2); if bb>nm, continue; end
    for q = 1:7
        rows = [rows; make_row(dataset,c,names(a)+" minus "+names(bb), ...
            metrics(q),values(a,q)-values(bb,q), ...
            boot(:,q,a)-boot(:,q,bb),nboot,seed)]; %#ok<AGROW>
    end
end
rows = [rows; make_row(dataset,c,"NR-TMM vs NR-full", ...
    "centroid_vector_error",hypot(values(1,1)-values(2,1), ...
    values(1,2)-values(2,2)),dc,nboot,seed)];
writetable(rows,file);
end

function rows = import_old(t)
rows = table();
for k=1:height(t)
    r=t(k,:);
    if startsWith(r.metric,"paired_landing"), continue; end
    if r.metric ~= "centroid_vector_error"
        rows=[rows; direct_row("table_iii",r.case_name,"NR-TMM",r.metric, ...
            r.model_a_value_nm,r.model_a_ci95_low_nm,r.model_a_ci95_high_nm, ...
            r.n_bootstrap,r.seed)]; %#ok<AGROW>
        rows=[rows; direct_row("table_iii",r.case_name,"NR-full",r.metric, ...
            r.model_b_value_nm,r.model_b_ci95_low_nm,r.model_b_ci95_high_nm, ...
            r.n_bootstrap,r.seed)]; %#ok<AGROW>
        label="NR-TMM minus NR-full";
    else
        label="NR-TMM vs NR-full";
    end
    rows=[rows; direct_row("table_iii",r.case_name,label,r.metric,r.estimate_nm, ...
        r.ci95_low_nm,r.ci95_high_nm,r.n_bootstrap,r.seed)]; %#ok<AGROW>
end
end

function r = make_row(ds,c,model,metric,value,boots,nboot,seed)
ci = quantile_linear(boots,[0.025 0.975]);
r = direct_row(ds,c,model,metric,value,ci(1),ci(2),nboot,seed);
end

function r = direct_row(ds,c,model,metric,value,lo,hi,nboot,seed)
r = table(string(ds),string(c),string(model),string(metric),value,lo,hi, ...
    nboot,seed,'VariableNames',{'dataset','case_name','model','metric', ...
    'value_nm','ci95_low_nm','ci95_high_nm','n_bootstrap','seed'});
end

function m = beam_metrics(xy)
c = mean(xy,1); r = hypot(xy(:,1)-c(1),xy(:,2)-c(2));
m = [c hypot(c(1),c(2)) sqrt(mean(r.^2)) ...
    2*quantile_linear(r,[.5 .9 .99])];
end

function q = quantile_linear(x,p)
x=sort(x(:)); pos=1+(numel(x)-1)*p; lo=floor(pos); hi=ceil(pos);
q=reshape(x(lo),size(p))+(pos-lo).*reshape(x(hi)-x(lo),size(p));
end

function f = load_comsol(file)
t=readmatrix(file,'FileType','text'); t=t(all(isfinite(t(:,[1 2 4])),2),:);
id=round(t(:,4)); assert(max(abs(t(:,4)-id))<1e-8);
[f.id,order]=sort(id); assert(numel(unique(f.id))==50000);
f.xy=t(order,1:2)*1e3;
end
