function run_t2_case06_axis_1um_20261002()
%RUN_T2_CASE06_AXIS_1UM_20261002 Independent near-axis Table III case.
    import com.comsol.model.util.*
    root = fileparts(fileparts(mfilename('fullpath')));
    out = fullfile(root,'RSI非空间电荷执行结果_20261001', ...
        '02_旧矩阵回归','Case06_axis_1um_20261002');
    if ~exist(out,'dir'), mkdir(out); end
    frozen = fullfile(root,'RSI非空间电荷执行结果_20261001', ...
        '02_旧矩阵回归','frozen_code');
    addpath(frozen,'-begin');
    source = fullfile(root,'01_intrinsic_error_initial.csv');
    [V,p] = load_comsol_initial_csv(source);
    assert(isequal(round(p.id(:)),(1:50000)'), ...
        'Frozen input pidx is incomplete.');
    [sys,err] = table_iii_tmm_config( ...
        '06_b_field_gradient_edge','legacy');
    sys.func_B_vec = @(X,Y,Z) [-3.*(X-1e-6); -3.*Y; ...
        sys.B0+6.*(Z-0.5.*sys.L)];
    assert(~isfield(sys,'b_gradient') || ~sys.b_gradient.enabled, ...
        'Unified B-gradient path must be disabled.');
    tmm = run_transfer_matrix_model(p,V,sys,err,numel(p.id));
    assert(~tmm.diagnostics.b_gradient.enabled, ...
        'Unified B-gradient operator was invoked.');
    tmm_xy = [tmm.X_final(:),tmm.Y_final(:)];
    writetable(table(p.id(:),tmm_xy(:,1),tmm_xy(:,2), ...
        'VariableNames',{'pidx','x_m','y_m'}), ...
        fullfile(out,'tmm_image_with_pidx.csv'));

    addpath('D:\Program Files\COMSOL\COMSOL64\Multiphysics\mli');
    user = getenv('COMSOL_REV_USER'); pass = getenv('COMSOL_REV_PASSWORD');
    assert(~isempty(user) && ~isempty(pass), ...
        'COMSOL LiveLink credentials are required.');
    mphstart('localhost',2036,user,pass);
    model_file = fullfile(root, ...
        '※磁场聚焦光电阴极模型（30纳米-毫米尺度）.mph');
    fprintf('Loading model: %s\n',model_file);
    model = mphload(model_file);
    cleanup = onCleanup(@() ModelUtil.remove(model.tag()));
    model.sol('sol1').clearSolutionData();
    cpt = model.component('comp1').physics('cpt');
    for tag = {'inl1','inl2','inl3','mf1'}
        cpt.feature(tag{1}).active(false);
    end
    cpt.feature('mf2').active(true);
    cpt.feature('ef1').active(true);
    cpt.feature('ef1').set('E',{'0','0','-1.267e7[V/m]'});
    cpt.feature('mf2').set('B',{'-3[T/m]*(x-1[um])', ...
        '-3[T/m]*y', ...
        '1.2[T]+6[T/m]*(z-0.001[um]-0.5[mm])'});
    cpt.prop('RelativisticCorrection').setIndex( ...
        'RelativisticCorrection',0,0);
    release = fullfile(root,'revision_comsol', ...
        'table_iii_nr_full_20260910', ...
        'canonical_intrinsic_particles_z_0p001um.txt');
    try
        cpt.feature.remove('rfd_T2_1um');
    catch
    end
    cpt.create('rfd_T2_1um','ReleaseFromDataFile');
    rfd = cpt.feature('rfd_T2_1um');
    rfd.set('Filename',release);
    rfd.set('icolp','0'); rfd.set('Rs','1');
    rfd.set('VelocitySpecification','SpecifyVelocity');
    rfd.set('InitialVelocity','FromFile');
    rfd.set('icolv','3'); rfd.set('rt',0);
    rfd.set('rc','7.125e-18[A]'); rfd.importData();
    model.study('std1').feature('time').set( ...
        'tlist','range(0,4e-15,4e-11)');
    solver = model.sol('sol1').feature('t1');
    solver.set('tlist','range(0,4e-15,4e-11)');
    solver.set('tstepsgenalpha','free');
    solver.set('maxstepconstraintgenalpha','const');
    solver.set('maxstepgenalpha','4e-15');
    write_settings(fullfile(out,'settings_before.txt'),model,model_file,release);
    diary(fullfile(out,'run.log'));
    started = datetime('now'); timer = tic;
    try
        model.sol('sol1').runAll();
        solve_seconds = toc(timer);
        exporter = model.result.export('data1');
        exporter.set('expr',{'cpt.pidx','qx','qy','qz','cpt.Ep', ...
            'cpt.vx','cpt.vy','cpt.vz'});
        exporter.set('unit',{'1','um','um','um','J', ...
            'm/s','m/s','m/s'});
        exporter.set('filename',fullfile(out,'initial_raw.csv'));
        exporter.run();
        section = model.result.dataset('cpl2');
        section.set('quickplane','xy');
        section.set('quickz','0.001000001[m]');
        poincare = model.result('pg2').feature('poma2');
        poincare.set('solnum','10001');
        poincare.set('looplevel',{'10001'});
        poincare.set('interp',{'4e-11'});
        poincare.set('usefixedpointsize',false);
        poincare.set('pointradiusexpr','cpt.pidx');
        exporter = model.result.export('plot1');
        exporter.set('filename',fullfile(out,'image_raw.csv'));
        exporter.run();
        write_settings(fullfile(out,'settings_after.txt'),model,model_file,release);
        fem = readmatrix(fullfile(out,'image_raw.csv'),'FileType','text');
        fem = fem(all(isfinite(fem(:,[1,2,4])),2),:);
        id = round(fem(:,4));
        assert(isequal(sort(id),(1:50000)') && ...
            all(abs(fem(:,4)-id)<1e-6),'Image pidx mismatch.');
        [~,order] = sort(id);
        fem_xy = fem(order,1:2)*1e-6;
        writetable(table(id(order),fem_xy(:,1),fem_xy(:,2), ...
            'VariableNames',{'pidx','x_m','y_m'}), ...
            fullfile(out,'fem_image_with_pidx.csv'));
        initial = readmatrix(fullfile(out,'initial_raw.csv'), ...
            'FileType','text');
        initial = initial(all(isfinite(initial(:,2:9)),2),:);
        assert(isequal(sort(round(initial(:,2))),(1:50000)'), ...
            'Initial pidx mismatch.');
        metrics = compare_metrics(tmm_xy,fem_xy);
        writetable(struct2table(metrics), ...
            fullfile(out,'case06_1um_metrics.csv'));
        fid = fopen(fullfile(out,'run_summary.txt'),'w');
        fprintf(fid,'status=success\nstarted=%s\nended=%s\n', ...
            char(started),char(datetime('now')));
        fprintf(fid,'solve_seconds=%.9g\nparticle_count=50000\n', ...
            solve_seconds);
        fclose(fid);
        fprintf('Case06 1um complete: centroid error %.6g nm, FW99 difference %+.6g nm.\n', ...
            metrics.centroid_error_nm,metrics.delta_FW99_nm);
    catch e
        fid = fopen(fullfile(out,'failure.txt'),'w');
        if fid >= 0
            fprintf(fid,'%s\n',getReport(e,'extended'));
            fclose(fid);
        end
        diary('off');
        rethrow(e);
    end
    diary('off');
    model.sol('sol1').clearSolutionData();
    clear cleanup;
end

function write_settings(path,model,model_file,release)
    cpt = model.component('comp1').physics('cpt');
    solver = model.sol('sol1').feature('t1');
    fid = fopen(path,'w');
    assert(fid >= 0,'Cannot write settings.');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid,'case=Case06_axis_1um\nmodel=%s\nrelease=%s\n', ...
        model_file,release);
    fprintf(fid,'version=%s\nhost=%s\n', ...
        char(com.comsol.model.util.ModelUtil.getComsolVersion()), ...
        getenv('COMPUTERNAME'));
    fprintf(fid,'E=%s\nB=%s\n', ...
        strjoin(string(cpt.feature('ef1').getStringArray('E')),';'), ...
        strjoin(string(cpt.feature('mf2').getStringArray('B')),';'));
    fprintf(fid,'relativistic_correction=%s\n', ...
        char(cpt.prop('RelativisticCorrection').getString( ...
        'RelativisticCorrection')));
    fprintf(fid,'max_step_s=%s\ntlist=%s\nrtol=%s\n', ...
        char(solver.getString('maxstepgenalpha')), ...
        char(solver.getString('tlist')),char(solver.getString('rtol')));
    fprintf(fid,'image_z=%s\n', ...
        char(model.result.dataset('cpl2').getString('quickz')));
    clear cleanup;
end

function m = compare_metrics(tmm,fem)
    a = beam_metrics(tmm); b = beam_metrics(fem);
    d = hypot(tmm(:,1)-fem(:,1),tmm(:,2)-fem(:,2))*1e9;
    m.case_name = "Case06_axis_1um";
    m.axis_offset_um = 1;
    m.N = size(tmm,1);
    m.tmm_cx_nm = a.cx; m.tmm_cy_nm = a.cy;
    m.fem_cx_nm = b.cx; m.fem_cy_nm = b.cy;
    m.centroid_error_nm = hypot(a.cx-b.cx,a.cy-b.cy);
    m.tmm_rms_nm = a.rms; m.fem_rms_nm = b.rms;
    m.rms_error_pct = 100*(a.rms-b.rms)/b.rms;
    m.tmm_FW99_nm = a.fw99; m.fem_FW99_nm = b.fw99;
    m.delta_FW99_nm = a.fw99-b.fw99;
    m.paired_rms_nm = sqrt(mean(d.^2));
    m.paired_p99_nm = quantile_linear(d,0.99);
    m.pass_centroid = m.centroid_error_nm < 1;
    m.pass_rms = abs(m.rms_error_pct) < 1;
    m.pass_FW99 = abs(m.delta_FW99_nm) < 1;
    m.pass_all = m.pass_centroid && m.pass_rms && m.pass_FW99;
end

function m = beam_metrics(xy)
    center = mean(xy,1);
    r = hypot(xy(:,1)-center(1),xy(:,2)-center(2));
    m.cx = center(1)*1e9; m.cy = center(2)*1e9;
    m.rms = sqrt(mean(r.^2))*1e9;
    m.fw99 = 2*quantile_linear(r,0.99)*1e9;
end

function q = quantile_linear(x,p)
    x = sort(x(:)); at = 1+(numel(x)-1)*p;
    lo = floor(at); hi = ceil(at);
    q = x(lo)+(at-lo)*(x(hi)-x(lo));
end
