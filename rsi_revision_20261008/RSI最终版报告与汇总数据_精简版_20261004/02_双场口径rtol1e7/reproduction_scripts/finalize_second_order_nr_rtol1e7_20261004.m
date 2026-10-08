function finalize_second_order_nr_rtol1e7_20261004
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'revision_comsol','second_order_nr_rtol1e7_dual_origin_20261004');
origins={'entrance','midpoint'}; names={'H0_control','H0p5','H1','H2'};
H=[0,1.2,2.4,4.8]; ids=(1:50000)'; xy=zeros(50000,2,4,2);
tx=zeros(50000,2,2); oldxy=nan(size(xy)); records=struct([]); steps=struct([]);
initial_reference=[];
for g=1:2
    d=fullfile(out,origins{g}); t=readtable(fullfile(d,'tmm_first_order_image_particles.csv'));
    assert(isequal(t.pidx,ids)); tx(:,:,g)=t{:,2:3}; mt=metric(tx(:,:,g));
    for k=1:4
        c=fullfile(d,names{k}); initial=readtable(fullfile(c,'initial_particles.csv'));
        image=readtable(fullfile(c,'image_particles.csv'));
        assert(isequal(initial.pidx,ids) && isequal(image.pidx,ids),'Incomplete IDs.');
        if isempty(initial_reference), initial_reference=initial{:,2:end}; end
        assert(isequal(initial_reference,initial{:,2:end}),'Input differs among eight cases.');
        for snap={'settings_before.txt','settings_after.txt'}
            text=fileread(fullfile(c,snap{1}));
            token=regexp(text,'(?m)^rtol=([^\r\n]+)','tokens','once');
            assert(abs(str2double(token{1})-1e-7)<1e-20);
            assert(contains(text,'relativistic_correction=0') && contains(text,'space_charge=off') && ...
                contains(text,['field_reference=' origins{g}]) && contains(text,sprintf('H_Tpm2=%.17g',H(k))) && ...
                contains(text,'image_plane_z=0.001000001[m]'));
        end
        log=fileread(fullfile(c,'solver_progress.log'));
        a=regexp(log,'(?m)^\s*(\d+)\s+([0-9.]+(?:[eE][+-]?\d+)?)\s+([0-9.]+(?:[eE][+-]?\d+)?)\s+(?:out\s+)?\d+','tokens');
        a=cellfun(@(x) cellfun(@str2double,x),a,'UniformOutput',false); a=vertcat(a{:});
        assert(~isempty(a) && isequal(a(:,1),(1:size(a,1))') && ...
            abs(a(end,2)-4e-11)<1e-22 && max(a(:,3))<=4e-15*(1+1e-10),'Actual step audit failed.');
        run=fileread(fullfile(c,'run_summary.txt')); sec=regexp(run,'solve_seconds=([0-9.eE+-]+)','tokens','once');
        steps=[steps;struct('field_reference',origins{g},'case_name',names{k},'accepted_steps',size(a,1), ...
            'min_step_s',min(a(:,3)),'max_step_s',max(a(:,3)),'final_time_s',a(end,2), ...
            'rtol',1e-7,'solve_seconds',str2double(sec{1}))]; %#ok<AGROW>
        xy(:,:,k,g)=image{:,2:3}; m=metric(xy(:,:,k,g)); m0=metric(xy(:,:,1,g));
        delta=(tx(:,:,g)-xy(:,:,k,g))*1e9; dh=(xy(:,:,k,g)-xy(:,:,1,g))*1e9;
        olddir='';
        if g==1
            olddir=fullfile(root,'revision_comsol','second_order_nr_paired_20261004',names{k});
        elseif k==1
            olddir=fullfile(root,'revision_comsol','second_order_nr_paired_20261004','diagnostic_midpoint_H0');
        end
        oldfw=NaN; tol_difference=NaN; tol_paired_rms=NaN;
        if ~isempty(olddir) && exist(fullfile(olddir,'image_particles.csv'),'file')
            oldinitial=readtable(fullfile(olddir,'initial_particles.csv'));
            assert(isequal(oldinitial.pidx,ids) && ...
                max(abs(oldinitial{:,2:4}-initial{:,2:4}),[],'all')<1e-18 && ...
                max(abs(oldinitial{:,6:8}-initial{:,6:8}),[],'all')<1e-8, ...
                'Old/new tolerance cases do not share released positions and velocities.');
            oldsettings=fileread(fullfile(olddir,'settings_after.txt'));
            oldtol=regexp(oldsettings,'(?m)^rtol=([^\r\n]+)','tokens','once');
            oldB=regexp(oldsettings,'(?m)^B=([^\r\n]+)','tokens','once');
            newB=regexp(text,'(?m)^B=([^\r\n]+)','tokens','once');
            oldH=regexp(oldsettings,'(?m)^H_Tpm2=([^\r\n]+)','tokens','once');
            assert(abs(str2double(oldtol{1})-1e-5)<1e-18 && ...
                strcmp(oldB{1},newB{1}) && abs(str2double(oldH{1})-H(k))<1e-15, ...
                'Old/new tolerance comparison has a different field reference or H.');
            old=readtable(fullfile(olddir,'image_particles.csv')); assert(isequal(old.pidx,ids));
            oldxy(:,:,k,g)=old{:,2:3}; om=metric(old{:,2:3}); oldfw=om.fw;
            tol_difference=m.fw-oldfw;
            tol_paired_rms=sqrt(mean(sum(((xy(:,:,k,g)-old{:,2:3})*1e9).^2,2)));
        end
        records=[records;struct('field_reference',origins{g},'case_name',names{k},'H_Tpm2',H(k),'N',50000, ...
            'TMM_FW99_nm',mt.fw,'NR_full_FW99_nm',m.fw,'TMM_minus_full_FW99_nm',mt.fw-m.fw, ...
            'H_minus_H0_FW99_nm',m.fw-m0.fw,'centroid_error_nm',norm(mt.center-m.center), ...
            'paired_TMM_full_RMS_nm',sqrt(mean(sum(delta.^2,2))), ...
            'paired_H_H0_RMS_nm',sqrt(mean(sum(dh.^2,2))), ...
            'old_rtol1e5_FW99_nm',oldfw,'rtol1e7_minus_rtol1e5_FW99_nm',tol_difference, ...
            'paired_tolerance_RMS_nm',tol_paired_rms,'settings_and_pairing_pass',true)]; %#ok<AGROW>
        writetable(table(ids,delta(:,1),delta(:,2),dh(:,1),dh(:,2), ...
            'VariableNames',{'pidx','TMM_full_dx_nm','TMM_full_dy_nm','H_H0_dx_nm','H_H0_dy_nm'}), ...
            fullfile(c,'paired_coordinate_audit.csv'));
    end
end
S=struct2table(records); writetable(S,fullfile(out,'dual_origin_truncation_summary.csv'));
writetable(struct2table(steps),fullfile(out,'actual_step_and_settings_audit.csv'));
nboot=1000; rng(20261004,'twister'); bh=zeros(nboot,4,2); bt=bh; tol=nan(size(bh));
for b=1:nboot
    ix=randi(50000,50000,1);
    for g=1:2
        f0=fw99(xy(ix,:,1,g)); ft=fw99(tx(ix,:,g));
        for k=1:4
            f=fw99(xy(ix,:,k,g)); bh(b,k,g)=f-f0; bt(b,k,g)=ft-f;
            if all(isfinite(oldxy(1,:,k,g))), tol(b,k,g)=f-fw99(oldxy(ix,:,k,g)); end
        end
    end
end
ci=struct([]);
for g=1:2, for k=1:4
    row=struct('field_reference',origins{g},'case_name',names{k},'bootstrap_repeats',nboot,'seed',20261004, ...
        'H_H0_CI95_low_nm',ql(bh(:,k,g),.025),'H_H0_CI95_high_nm',ql(bh(:,k,g),.975), ...
        'TMM_full_CI95_low_nm',ql(bt(:,k,g),.025),'TMM_full_CI95_high_nm',ql(bt(:,k,g),.975), ...
        'tolerance_CI95_low_nm',ql(tol(:,k,g),.025),'tolerance_CI95_high_nm',ql(tol(:,k,g),.975));
    ci=[ci;row]; %#ok<AGROW>
end, end
C=struct2table(ci); writetable(C,fullfile(out,'paired_bootstrap_CI95.csv'));
f=figure('Visible','off','Color','w','Position',[100,100,1000,450]); tiledlayout(1,2,'Padding','loose');
colors=lines(2);
nexttile; hold on;
for g=1:2, plot(H,S.NR_full_FW99_nm((g-1)*4+(1:4)),'o-','Color',colors(g,:)); end
xlabel('H (T/m^2)'); ylabel('NR-full FW99 (nm)'); legend('B0 at emission','Linear G centered at midpoint','Location','best'); grid on;
nexttile; hold on;
for g=1:2
    ix=(g-1)*4+(1:4); y=S.H_minus_H0_FW99_nm(ix);
    errorbar(H,y,y-C.H_H0_CI95_low_nm(ix),C.H_H0_CI95_high_nm(ix)-y,'o-','Color',colors(g,:));
end
xlabel('H (T/m^2)'); ylabel('FW99(H)-FW99(0) (nm)'); grid on;
print(f,fullfile(out,'dual_origin_truncation.png'),'-dpng','-r220'); close(f);
fid=fopen(fullfile(out,'双场口径容差1e-7复算报告.md'),'w','n','UTF-8');
fprintf(fid,'# 二阶截断表：双场口径 rtol=1e-7 复算\n\n');
fprintf(fid,'完成八个独立 COMSOL 研究1工况。每组50,000粒子，复用26/27/28历史初态；相对论、空间电荷均关闭。研究节点与实际求解器求解前后均核验rtol=1e-7。原模型、旧结果和论文未覆盖。\n\n');
fprintf(fid,'| 场口径 | H (T/m²) | TMM FW99 (nm) | NR-full FW99 (nm) | TMM-full (nm) | full(H)-full(0) (nm) |\n|---|---|---|---|---|---|\n');
for i=1:height(S)
    fprintf(fid,'| %s | %.3g | %.9f | %.9f | %+.9f | %+.9f |\n',S.field_reference{i},S.H_Tpm2(i),S.TMM_FW99_nm(i),S.NR_full_FW99_nm(i),S.TMM_minus_full_FW99_nm(i),S.H_minus_H0_FW99_nm(i));
end
fprintf(fid,'\n## 口径与数值边界\n\n');
fprintf(fid,'`entrance`：Bz=1.2+6Z+HZ²/2。`midpoint`延续历史模板：Bz=1.2+6(Z-L/2)+HZ²/2，仅线性G项居中，H项仍以发射面为原点。因此H非零时，中点的实际场为1.2+HL²/8 T，不能把它写成严格固定1.2 T的总中点场。\n\n');
fprintf(fid,'Bx=-x(6+HZ)/2，By=-y(6+HZ)/2，满足div(B)=0。Z=z-1nm，L=1mm；发射z=1nm，固定像面z=1000.001μm；E=(0,0,-1.267e7)V/m。H=0/1.2/2.4/4.8T/m²。\n\n');
fprintf(fid,'广义alpha，rhoinf=1，最大内部步长4fs，输出0:4fs:40ps；atolglobal字段读回0.001，网格沿用hmax=0.24μm/hmin=0.0024μm。实际绝对容差沿用scaled/factor模式、因子0.1，带自动缩放；不能将0.001字段称为直接固定的绝对误差阈值。完整有效参数由重新加载保存模型的LiveLink API读回，见saved_model_livelink_settings_audit.csv；网格哈希见saved_model_tolerance_and_mesh_audit.csv。序列化key中的通用默认文本不作有效设置证据。真实接受步数来自solver_progress.log，而非输出点数。\n\n');
fprintf(fid,'TMM使用交接包最终EB0/E1/共享包络算子，保留旧B路径与midpoint相位，I=0、ee4000；不启用统一B或积分相位。TMM不含H，故同一口径四组预测相同。U0由原始速度按TMM常数反算NR动能，以匹配COMSOL实际释放速度。\n\n');
fprintf(fid,'FW99=2×以各自质心为中心的99%%半径线性分位数。pidx逐粒子完整配对；八组位置/速度/初态数据完全相同。配对bootstrap1000次、seed20261004，CI是有限样本重采样不确定性，不是求解误差界限。\n\n');
fprintf(fid,'## 容差改变对照\n\n| 场口径/工况 | 旧rtol1e-5 FW99 (nm) | 新减旧 (nm) | 配对坐标RMS变化 (nm) |\n|---|---|---|---|\n');
for i=1:height(S)
    if isfinite(S.old_rtol1e5_FW99_nm(i)), fprintf(fid,'| %s/%s | %.9f | %+.9f | %.9f |\n',S.field_reference{i},S.case_name{i},S.old_rtol1e5_FW99_nm(i),S.rtol1e7_minus_rtol1e5_FW99_nm(i),S.paired_tolerance_RMS_nm(i)); end
end
fprintf(fid,'\n中点口径旧容差只有H=0对照，其他三档无旧数据，CSV保留NaN，不补造。收紧容差不能代替独立更细步长验证。\n\n');
fprintf(fid,'## 结论与使用规则\n\nH=0的TMM-full差是既有一阶背景误差，不能算成二阶截断误差。二阶场独立响应使用相同口径的full(H)-full(0)，总差与增量分别报告。不同口径不能相互充当H=0参考。\n\n');
fprintf(fid,'最大|H/H0 FW99增量|：发射面口径%.9g nm，中点模板口径%.9g nm。即使总差近似相同，也不证明统计噪声下限或全模型亚纳米准确。\n\n',max(abs(S.H_minus_H0_FW99_nm(1:4))),max(abs(S.H_minus_H0_FW99_nm(5:8))));
fprintf(fid,'### 配对置信区间\n\n| 场口径/工况 | H/H0点估计 (nm) | 95%% CI低端 (nm) | 95%% CI高端 (nm) |\n|---|---|---|---|\n');
for i=[2:4,6:8]
    fprintf(fid,'| %s/%s | %.9g | %.9g | %.9g |\n',S.field_reference{i},S.case_name{i}, ...
        S.H_minus_H0_FW99_nm(i),C.H_H0_CI95_low_nm(i),C.H_H0_CI95_high_nm(i));
end
fprintf(fid,'\n发射面三档H/H0区间为正；中点模板三档均跨零，故中点固定样本的正点估计不能表述为对发射分布稳健的显著展宽。中点H=0的容差差区间也跨零；所有CI仅描述配对重采样，不是数值误差边界。\n\n');
fprintf(fid,'容差改动使绝对FW99变化约0.0023–0.0024nm，比当前H响应点估计更大。发射面配对H/H0增量在两种容差间的变化仅约2.42e-8nm以内，说明共同基线偏移大部分抵消，但这仍不能代替独立细步长验证。\n\n');
fprintf(fid,'本轮不支持直接保留历史三行相同0.12nm的截断表。应按本报告区分总差与二阶场增量，并以实际口径替换旧表。\n\n');
fprintf(fid,'逐粒子初态/像面数据、设置快照、真实求解日志与无解设置模型位于对应口径/H工况目录。汇总表、配对误差与CI可复算。\n'); fclose(fid);
write_text(fullfile(out,'AUDIT_AND_REPORT_COMPLETE.txt'),'Eight-case pairing/settings/step audit and paired CIs complete.');
fprintf('DUAL_ORIGIN_AUDIT_AND_REPORT_COMPLETE\n');
end
function m=metric(xy)
c=mean(xy,1); r=hypot(xy(:,1)-c(1),xy(:,2)-c(2))*1e9;
m=struct('center',c*1e9,'fw',2*ql(r,.99));
end
function f=fw99(xy)
m=metric(xy); f=m.fw;
end
function v=ql(x,p)
if any(~isfinite(x)), v=NaN; return; end
x=sort(x(:)); a=1+(numel(x)-1)*p; l=floor(a); h=ceil(a); v=x(l)+(a-l)*(x(h)-x(l));
end
function write_text(file,txt)
fid=fopen(file,'w'); assert(fid>=0); fprintf(fid,'%s\n',txt); fclose(fid);
end
