function plot_reviewer7_sensitivity()
% Plot published statistics without rerunning COMSOL or changing physics.
here = fileparts(mfilename('fullpath')); root = fileparts(here);
data = fullfile(root,'RSI仿真结果总整理_20261002');
set(groot,'defaultAxesFontName','Times New Roman', ...
    'defaultTextFontName','Times New Roman');
colors = [0.12 0.32 0.67; 0.80 0.25 0.12];
s=readtable(fullfile(data,'04_统计与图件','mc_convergence_summary.csv'));
seeds=readtable(fullfile(data,'04_统计与图件','mc_seed_summary.csv'));
metrics={'FW90','FW99'};
fig=figure('Visible','off','Color','w','Units','inches', ...
    'Position',[.5 .5 7.2 4.6]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for j=1:2
    ax=nexttile(tl,j); hold(ax,'on');
    meanv=s.([metrics{j} '_mean_nm']);
    lo=s.([metrics{j} '_ci95_low_nm']); hi=s.([metrics{j} '_ci95_high_nm']);
    for k=1:height(s)
        v=seeds.(sprintf('%s_nm',metrics{j}))(seeds.N==s.N(k));
        x=s.N(k)*10.^linspace(-.035,.035,numel(v))';
        scatter(ax,x,v,11,colors(j,:),'filled','MarkerFaceAlpha',.30);
    end
    errorbar(ax,s.N,meanv,meanv-lo,hi-meanv,'o-', ...
        'Color',colors(j,:),'MarkerFaceColor',colors(j,:), ...
        'LineWidth',1.1,'MarkerSize',4,'CapSize',5);
    set(ax,'XScale','log'); xlim(ax,[650 8e7]);
    allv = seeds.(sprintf('%s_nm',metrics{j}));
    padding = .08*(max(allv)-min(allv));
    ylim(ax,[min([allv;lo])-padding max([allv;hi])+padding]);
    ylabel(ax,[metrics{j} ' (nm)']);
    title(ax,sprintf('(%c) %s; ten seeds per count','a'+j-1,metrics{j}));
    xticks(ax,s.N); ax.XTickLabel=[]; grid(ax,'on');
    ax=nexttile(tl,j+2); hold(ax,'on');
    ref=meanv(end);
    patch(ax,[650 8e7 8e7 650], ...
        [lo(end)-ref lo(end)-ref hi(end)-ref hi(end)-ref], ...
        [.88 .88 .88],'EdgeColor','none');
    yline(ax,0,':','Color',[.35 .35 .35]);
    errorbar(ax,s.N,meanv-ref,meanv-lo,hi-meanv,'o-', ...
        'Color',colors(j,:),'MarkerFaceColor',colors(j,:), ...
        'LineWidth',1.1,'MarkerSize',4,'CapSize',5);
    set(ax,'XScale','log'); xlim(ax,[650 8e7]); xticks(ax,s.N);
    xticklabels(ax,{'10^3','10^4','10^5','10^6','10^7','5x10^7'});
    xlabel(ax,'Particle count'); ylabel(ax,['Mean ' metrics{j} ' - reference (nm)']);
    title(ax,sprintf('(%c) Relative to largest-count mean','c'+j-1)); grid(ax,'on');
end
save_figure(fig,root,'fig_particle_count_reviewer7');

s=readtable(fullfile(data,'03_时间步与网格','mesh_M1_M5_case_summary.csv'), ...
    'TextType','string');
ci=readtable(fullfile(data,'03_时间步与网格','mesh_M1_M5_bootstrap_ci.csv'), ...
    'TextType','string');
inverse_hmax = 1 ./ s.hmax_um;
assert(all(diff(inverse_hmax)>0));
mesh_xlim = [inverse_hmax(1)/1.25 inverse_hmax(end)*1.25];
mesh_ticklabels = arrayfun(@(x) sprintf('%.3f', x), inverse_hmax, ...
    'UniformOutput', false);
fig=figure('Visible','off','Color','w','Units','inches', ...
    'Position',[.5 .5 3.5 4.45]);
tl=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); hold(ax,'on');
for j=1:2
    v=s.(['delta_' metrics{j} '_nm']); lo=zeros(5,1); hi=lo;
    for k=1:5
        r=ci(ci.case_id==s.case_id(k) & ci.metric=="delta_"+metrics{j}+"_nm",:);
        if ~isempty(r), lo(k)=r.ci95_lower_nm; hi(k)=r.ci95_upper_nm; end
    end
    markers = {'o','s'}; sizes = [5 3.5];
    face = 'w'; if j==2, face=colors(j,:); end
    errorbar(ax,inverse_hmax,v,v-lo,hi-v,[markers{j} '-'], ...
        'Color',colors(j,:),'MarkerFaceColor',face, ...
        'MarkerSize',sizes(j),'LineWidth',1.1,'DisplayName',metrics{j});
end
yline(ax,0,':','HandleVisibility','off'); ylim(ax,[-.01 .01]);
xlim(ax,mesh_xlim); set(ax,'XScale','log');
xticks(ax,inverse_hmax); xticklabels(ax,mesh_ticklabels);
ylabel(ax,'Width - M5 width (nm)'); title(ax,'(a) Analytic-field mesh sensitivity');
legend(ax,'Location','northeast','Box','off','FontSize',8);
text(ax,sqrt(prod(mesh_xlim)), -.006,'No resolved mesh dependence', ...
    'HorizontalAlignment','center','FontSize',8); grid(ax,'on');
set(ax,'XMinorGrid','off');
ax=nexttile(tl); hold(ax,'on');
loglog(ax,inverse_hmax,s.total_elements,'o-','Color',[.25 .45 .35], ...
    'MarkerFaceColor',[.25 .45 .35],'LineWidth',1.1,'MarkerSize',4);
set(ax,'XScale','log','YScale','log'); ylim(ax,[1e6 4e8]); xlim(ax,mesh_xlim);
label_scale = [.72 1.5 1.9 1.25 1.25];
for k=1:5
    align='center'; if k>=3, align='right'; end
    if k==1, align='left'; end
    text(ax,inverse_hmax(k),s.total_elements(k)*label_scale(k), ...
        sprintf('M%d: %.2f M',k,s.total_elements(k)/1e6), ...
        'HorizontalAlignment',align,'FontSize',7.5);
end
xticks(ax,inverse_hmax); xticklabels(ax,mesh_ticklabels);
xlabel(ax,{'Reciprocal nominal maximum size', ...
    ['1/h_{max} (' char(956) 'm^{-1})']});
ylabel(ax,'Actual total elements');
title(ax,'(b) Realized mesh sizes'); grid(ax,'on');
set(ax,'XMinorGrid','off','YMinorGrid','off');
save_figure(fig,root,'fig_mesh_sensitivity_reviewer7');

s=readtable(fullfile(data,'03_时间步与网格','timestep_NR_metrics.csv'));
ci=readtable(fullfile(data,'03_时间步与网格','timestep_NR_bootstrap_ci.csv'), ...
    'TextType','string');
fig=figure('Visible','off','Color','w','Units','inches', ...
    'Position',[.5 .5 7.2 2.9]);
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
for panel=1:2
    ax=nexttile(tl); hold(ax,'on');
    yline(ax,0,':','HandleVisibility','off','Color',[.35 .35 .35]);
    for j=1:2
        v=s.(['delta_' metrics{j} '_nm']); lo=zeros(height(s),1); hi=lo;
        for k=1:height(s)
            r=ci(ci.equivalent_steps==s.equivalent_steps(k) & ...
                ci.metric=="delta_"+metrics{j}+"_nm",:);
            if s.equivalent_steps(k)==80000
                assert(isempty(r)); lo(k)=0; hi(k)=0; continue;
            end
            assert(height(r)==1);
            lo(k)=r.ci_low_nm; hi(k)=r.ci_high_nm;
        end
        errorbar(ax,s.equivalent_steps,v,v-lo,hi-v,'o-', ...
            'Color',colors(j,:),'MarkerFaceColor',colors(j,:), ...
            'LineWidth',1.1,'MarkerSize',4,'CapSize',5,'DisplayName',metrics{j});
    end
    set(ax,'XScale','log'); xlabel(ax,'Equivalent maximum-step count');
    ylabel(ax,'Width - 80k width (nm)'); grid(ax,'on');
    if panel==1
        xlim(ax,[80 1e5]); ylim(ax,[-1.15 .22]);
        title(ax,'(a) Full scan; paired 95% intervals');
        legend(ax,'Location','southeast','Box','off','FontSize',8);
    else
        xlim(ax,[800 1e5]); ylim(ax,[-.07 .21]);
        xticks(ax,[1000 2000 5000 10000 20000 40000 80000]);
        xticklabels(ax,{'1k','2k','5k','10k','20k','40k','80k'});
        title(ax,'(b) Fine-step region (non-monotonic)');
    end
end
save_figure(fig,root,'fig_timestep_sensitivity_reviewer7');
end

function save_figure(fig,root,name)
set(findall(fig,'-property','FontName'),'FontName','Times New Roman');
set(findall(fig,'Type','axes'),'FontSize',9,'LineWidth',.65,'Box','on');
file=fullfile(root,'figures',name);
exportgraphics(fig,[file '.pdf'],'ContentType','vector');
exportgraphics(fig,[file '.png'],'Resolution',300);
close(fig); fprintf('Saved %s\n',file);
end
