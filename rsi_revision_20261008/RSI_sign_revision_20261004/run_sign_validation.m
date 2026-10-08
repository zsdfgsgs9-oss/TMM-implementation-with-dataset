function run_sign_validation(nboot)
if nargin<1, nboot=5000; end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
runtime=fullfile(here,'runtime'); oldbr=fullfile(here,'old_br');
addpath(runtime,'-begin');
clear run_transfer_matrix_model calc_M_EB0 calc_M_E1 calc_M_Br1 table_iii_tmm_config
assert(strcmp(which('run_transfer_matrix_model'),fullfile(runtime,'run_transfer_matrix_model.m')));
base=fullfile(root,'RSI大修仿真统一归档_20260911','03_COMSOL完整结果', ...
    'revision_comsol','table_iii_nr_full_20260910');
initial=fullfile(root,'RSI_revision_simulation_package_20260901', ...
    '07_source_inputs','intrinsic_pair','01_intrinsic_error_initial.csv');
[v,p]=load_comsol_initial_csv(initial); p.epsilon_ee=zeros(size(p.U0));
names=["01_intrinsic_error","02_em_field_angle_error","03_e_field_zero_order", ...
    "04_b_field_zero_order","05_b_field_gradient","06_b_field_gradient_axis_1um", ...
    "07_e_field_gradient"];
stages=["before","electric_only","charge_only","both"];
ablations=table(); final=table();
intervals=readtable(fullfile(root,'RSI_reviewer7_20261004','table_iii_intervals.csv'), ...
    'TextType','string');
for k=1:7
    pk=p; vk=v;
    if k==6
        case6=fullfile(root,'RSI非空间电荷执行结果_20261001','02_旧矩阵回归', ...
            'Case06_axis_1um_20261002');
        [vk,pk]=load_comsol_initial_csv(fullfile(case6,'initial_raw.csv'));
        pk.epsilon_ee=zeros(size(pk.U0));
        f=sortrows(readtable(fullfile(case6,'fem_image_with_pidx.csv')),'pidx');
        fid=f.pidx; fullxy=[f.x_m f.y_m]*1e9;
    else
        reference=fullfile(base,char(names(k)),'rel_off','poincare_design_plane_with_pidx.csv');
        if k==1
            reference=fullfile(fileparts(base),'baseline_three_model_20260904', ...
                'NR_full','poincare_zrel_1000um_with_pidx.csv');
        end
        raw=readmatrix(reference,'FileType','text');
        raw=raw(all(isfinite(raw(:,[1 2 4])),2),:);
        [fid,ix]=sort(round(raw(:,4))); fullxy=raw(ix,1:2)*1e3;
    end
    [ids,ix]=sort(pk.id(:)); assert(isequal(ids,fid) && numel(unique(ids))==50000);
    allxy=cell(4,1);
    for s=1:4
        [sys,err]=table_iii_tmm_config(names(k));
        if s==1 || s==3
            E=sys.U/sys.L; ev=sys.func_E_vec;
            sys.func_E_vec=@(X,Y) flip_axial(ev(X,Y),E);
        end
        if s<3
            addpath(oldbr,'-begin'); clear calc_M_Br1
        else
            if contains(path,oldbr), rmpath(oldbr); end
            clear calc_M_Br1
        end
        out=run_transfer_matrix_model(pk,vk,sys,err,50000);
        assert(~out.diagnostics.b_gradient.enabled);
        allxy{s}=[out.X_final(ix).' out.Y_final(ix).']*1e9;
        mt=metrics(allxy{s}); mf=metrics(fullxy);
        ablations=[ablations;table(names(k),stages(s), ...
            norm(mt(1:2)-mf(1:2)),mt(4),mt(7),mt(7)-mf(7), ...
            'VariableNames',{'case_name','stage','centroid_error_nm','RMS_nm','FW99_nm','delta_FW99_nm'})]; %#ok<AGROW>
    end
    xy=allxy{4}; mt=metrics(xy); mf=metrics(fullxy);
    if k==6
        previous=sortrows(readtable(fullfile(case6,'tmm_image_with_pidx.csv')),'pidx');
        assert(isequal(previous.pidx,ids));
        assert(max(hypot(allxy{1}(:,1)-previous.x_m*1e9, ...
            allxy{1}(:,2)-previous.y_m*1e9))<2e-5,'Case 06 old replay failed.');
    end
    final=[final;table(names(k),norm(mt(1:2)-mf(1:2)),mt(4),mf(4), ...
        mt(5),mf(5),mt(6),mf(6),mt(7),mf(7),mt(7)-mf(7), ...
        'VariableNames',{'case_name','centroid_error_nm','TMM_RMS_nm','full_RMS_nm', ...
        'TMM_FW50_nm','full_FW50_nm','TMM_FW90_nm','full_FW90_nm', ...
        'TMM_FW99_nm','full_FW99_nm','delta_FW99_nm'})]; %#ok<AGROW>
    writetable(table(ids,xy(:,1),xy(:,2),fullxy(:,1),fullxy(:,2), ...
        'VariableNames',{'pidx','tmm_x_nm','tmm_y_nm','nr_full_x_nm','nr_full_y_nm'}), ...
        fullfile(here,sprintf('case%02d_coordinates.csv',k)));
    if k==5 || k==6
        if k==5, assert(abs(metrics(allxy{1})*[zeros(6,1);1]-39.0541169682701)<1e-6); end
        seed=2026100402+(k==6); boot=zeros(nboot,7,2); dc=zeros(nboot,1);
        for b=1:nboot
            stream=RandStream('Threefry','Seed',seed); stream.Substream=b;
            index=randi(stream,50000,50000,1);
            boot(b,:,1)=metrics(xy(index,:)); boot(b,:,2)=metrics(fullxy(index,:));
            dc(b)=norm(boot(b,1:2,1)-boot(b,1:2,2));
            if mod(b,1000)==0, fprintf('%s bootstrap %d/%d\n',names(k),b,nboot); end
        end
        one=table(); labels=["cx","cy","centroid_magnitude","RMS_radius","FW50","FW90","FW99"];
        model=["NR-TMM","NR-full","NR-TMM minus NR-full"];
        for j=1:3
            for q=1:7
                a=mt(q)-mf(q); bs=boot(:,q,1)-boot(:,q,2);
                if j==1, a=mt(q); bs=boot(:,q,1); end
                if j==2, a=mf(q); bs=boot(:,q,2); end
                one=[one; interval(names(k),model(j),labels(q),a,bs,nboot,seed)]; %#ok<AGROW>
            end
        end
        one=[one;interval(names(k),"NR-TMM vs NR-full","centroid_vector_error", ...
            norm(mt(1:2)-mf(1:2)),dc,nboot,seed)];
        intervals=[intervals(intervals.case_name~=names(k),:);one];
        writetable(one,fullfile(here,sprintf('case%02d_intervals.csv',k)));
    else
        old=intervals(intervals.case_name==names(k) & intervals.model=="NR-TMM" & intervals.metric=="FW99",:);
        assert(abs(old.value_nm-mt(7))<1e-7,'Unaffected case changed.');
        assert(max(abs(allxy{1}-xy),[],'all')<1e-9);
    end
    writetable(final,fullfile(here,'table_iii_metrics.csv'));
    writetable(ablations,fullfile(here,'sign_ablations.csv'));
    writetable(sortrows(intervals,{'case_name','metric','model'}),fullfile(here,'table_iii_intervals.csv'));
    disp(final(end,:));
end
run_sign_stripes(nboot);
end

function ev=flip_axial(ev,E)
ev(3,:)=E;
end
function m=metrics(xy)
c=mean(xy,1); r=hypot(xy(:,1)-c(1),xy(:,2)-c(2));
m=[c,norm(c),sqrt(mean(r.^2)),2*qlinear(r,[.5 .9 .99])];
end
function q=qlinear(x,p)
x=sort(x(:)); pos=1+(numel(x)-1)*p;
q=reshape(x(floor(pos)),size(p))+(pos-floor(pos)).*reshape(x(ceil(pos))-x(floor(pos)),size(p));
end
function r=interval(c,model,metric,a,bs,nboot,seed)
ci=qlinear(bs,[.025 .975]);
r=table("table_iii",c,model,metric,a,ci(1),ci(2),nboot,seed, ...
    'VariableNames',{'dataset','case_name','model','metric','value_nm', ...
    'ci95_low_nm','ci95_high_nm','n_bootstrap','seed'});
end
