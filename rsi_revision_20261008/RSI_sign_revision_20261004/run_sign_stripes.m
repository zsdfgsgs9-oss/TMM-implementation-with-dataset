function run_sign_stripes(nboot)
% Replay the signed published operators against existing COMSOL coordinates.
if nargin<1, nboot=5000; end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
source=fullfile(root,'Fig10最小工况重跑结果_20261002');
runtime=fullfile(here,'runtime'); addpath(runtime,'-begin');
clear calc_M_EB0 calc_M_E1 calc_M_Br1 calc_common_physics run_transfer_matrix_model
assert(strcmp(which('run_transfer_matrix_model'),fullfile(runtime,'run_transfer_matrix_model.m')));
sys=struct('B0',1.2,'L',1e-3,'U',12.67e3,'N_turn',1, ...
    'physics_model','nonrelativistic');
E=sys.U/sys.L; G=6;
sys.func_E_vec=@(X,Y)[.8e-3*E*ones(size(X));.6e-3*E*ones(size(Y));-E*ones(size(X))];
sys.func_B_vec=@(X,Y,Z)[-.5*G*X;-.5*G*Y;sys.B0+G*(Z-.5*sys.L)];
err=struct('G',G,'delta_E',1.267e4,'delta_B',1.2e-3,'z1',0, ...
    'z2',sys.L,'I_beam',0,'R0',18e-9,'ee_method','gauss');
err.func_dEz=@(x,y,z)6.335e7*z.*ones(size(x));
names={'11_center','13_axis_minus_1um'};
inputs={'11_system_stripe_center','13_system_stripe_1um'};
offsets=[0 1e-6];
previous=readtable(fullfile(root,'RSI_reviewer8_20261004','stripe_paired_coordinates.csv'), ...
    'TextType','string');
xy=cell(2,2); v=cell(2,1); p=cell(2,1); old=cell(2,1);
for k=1:2
    raw=readmatrix(fullfile(source,'inputs',[inputs{k} '_initial.csv']), ...
        'FileType','text');
    raw=raw(all(isfinite(raw(:,1:7)),2),1:7);
    assert(size(raw,1)==50000 && isequal(round(raw(:,1)),(1:50000)'));
    p{k}=struct('id',round(raw(:,1)).','U0',raw(:,4).'/1.602e-19, ...
        'alpha',atan2(hypot(raw(:,5),raw(:,6)),raw(:,7)).', ...
        'phi',atan2(raw(:,6),raw(:,5)).','epsilon_ee',zeros(1,50000));
    v{k}=zeros(3,1,50000);
    v{k}(1,1,:)=raw(:,2).'*1e-6+offsets(k);
    v{k}(2,1,:)=raw(:,3).'*1e-6; v{k}(3,1,:)=1;
    result=run_transfer_matrix_model(p{k},v{k},sys,err,50000);
    assert(~result.diagnostics.b_gradient.enabled);
    old{k}=[result.X_final(:)-offsets(k),result.Y_final(:)]*1e9;
    pair=sortrows(previous(previous.case_name==string(names{k}),:),'pidx');
    assert(isequal(pair.pidx,(1:50000)'));
    old{k}=[pair.tmm_x_nm pair.tmm_y_nm];
    f=sortrows(readtable(fullfile(source,names{k},'image_particles.csv')),'pidx');
    assert(isequal(f.pidx,pair.pidx));
    xy{k,2}=[f.x_m f.y_m]*1e9;
    assert(max(abs(xy{k,2}-[pair.nr_full_x_nm pair.nr_full_y_nm]),[],'all')<1e-9);
end
initial=sortrows(readtable(fullfile(source,names{1},'initial_particles.csv')),'pidx');
other=sortrows(readtable(fullfile(source,names{2},'initial_particles.csv')),'pidx');
assert(isequaln(initial,other));
centers=(-2:2)*72e-9;
[distance,labels]=min(abs(initial.x_m-centers),[],2);
assert(max(distance)<18e-9+1e-14 && all(accumarray(labels,1)==10000));
coordinates=table(); estimates=zeros(2,2,8);
for k=1:2
    result=run_transfer_matrix_model(p{k},v{k},sys,err,50000);
    assert(~result.diagnostics.b_gradient.enabled);
    xy{k,1}=[result.X_final(:)-offsets(k),result.Y_final(:)]*1e9;
    d=hypot(xy{k,1}(:,1)-xy{k,2}(:,1),xy{k,1}(:,2)-xy{k,2}(:,2));
    coordinates=[coordinates;table(repmat(string(names{k}),50000,1), ...
        initial.pidx,labels,xy{k,1}(:,1),xy{k,1}(:,2), ...
        xy{k,2}(:,1),xy{k,2}(:,2),d,'VariableNames', ...
        {'case_name','pidx','stripe_id','tmm_x_nm','tmm_y_nm', ...
        'nr_full_x_nm','nr_full_y_nm','paired_error_nm'})]; %#ok<AGROW>
    for j=1:2, estimates(k,j,:)=stripe_metrics(xy{k,j},labels); end
    fprintf('%s maximum sign-correction coordinate update %.9f nm\n',names{k}, ...
        max(hypot(xy{k,1}(:,1)-old{k}(:,1),xy{k,1}(:,2)-old{k}(:,2))));
end
writetable(coordinates,fullfile(here,'stripe_paired_coordinates.csv'));
metric_names=["cx_nm","cy_nm","RMS_radius_nm","FW99_nm", ...
    "pitch_nm","gaussian_equivalent_width_nm","stripe_FW90_nm","rotation_deg"];
boot=zeros(nboot,2,2,8); error_boot=zeros(nboot,2,2);
seed=2026100408; checkpoint=fullfile(here,'stripe_bootstrap_eb0_e1_checkpoint.mat');
done=0;
if isfile(checkpoint)
    saved=load(checkpoint);
    assert(saved.nboot==nboot && saved.seed==seed && isequal(saved.estimates,estimates));
    boot=saved.boot; error_boot=saved.error_boot; done=saved.done;
end
for b=done+1:nboot
    stream=RandStream('Threefry','Seed',seed); stream.Substream=b;
    index=randi(stream,50000,50000,1); sampled_labels=labels(index);
    for k=1:2
        for j=1:2
            boot(b,k,j,:)=stripe_metrics(xy{k,j}(index,:),sampled_labels);
        end
        error_boot(b,k,1)=hypot(boot(b,k,1,1)-boot(b,k,2,1), ...
            boot(b,k,1,2)-boot(b,k,2,2));
        d=xy{k,1}(index,:)-xy{k,2}(index,:);
        error_boot(b,k,2)=sqrt(mean(sum(d.^2,2)));
    end
    if mod(b,1000)==0 || b==nboot
        done=b; save(checkpoint,'boot','error_boot','done','seed','nboot','estimates');
        fprintf('Joint paired stripe bootstrap: %d/%d\n',b,nboot);
    end
end
metrics=table(); intervals=table(); response=table();
for k=1:2
    for q=1:8
        a=estimates(k,1,q); f=estimates(k,2,q);
        metrics=[metrics;table(string(names{k}),metric_names(q),a,f,a-f, ...
            'VariableNames',{'case_name','metric','NR_TMM','NR_full','difference'})]; %#ok<AGROW>
        for j=1:3
            model="NR-TMM"; value=a; bs=boot(:,k,1,q);
            if j==2, model="NR-full"; value=f; bs=boot(:,k,2,q); end
            if j==3
                model="NR-TMM minus NR-full"; value=a-f;
                bs=boot(:,k,1,q)-boot(:,k,2,q);
            end
            intervals=[intervals;interval_row(names{k},metric_names(q), ...
                model,value,bs,nboot,seed)]; %#ok<AGROW>
        end
    end
    d=xy{k,1}-xy{k,2}; distances=hypot(d(:,1),d(:,2));
    error_values=[norm(mean(d,1)),sqrt(mean(distances.^2))];
    for q=1:2
        mn=["centroid_vector_error_nm","paired_landing_RMS_nm"];
        metrics=[metrics;table(string(names{k}),mn(q),nan,nan,error_values(q), ...
            'VariableNames',{'case_name','metric','NR_TMM','NR_full','difference'})]; %#ok<AGROW>
        intervals=[intervals;interval_row(names{k},mn(q),"NR-TMM vs NR-full", ...
            error_values(q),error_boot(:,k,q),nboot,seed)]; %#ok<AGROW>
    end
end
for q=1:2
    metric=["axis_response_x_nm","axis_response_y_nm"];
    for j=1:3
        t=estimates(2,1,q)-estimates(1,1,q);
        f=estimates(2,2,q)-estimates(1,2,q);
        model="NR-TMM"; value=t; bs=boot(:,2,1,q)-boot(:,1,1,q);
        if j==2, model="NR-full"; value=f; bs=boot(:,2,2,q)-boot(:,1,2,q); end
        if j==3
            model="NR-TMM minus NR-full"; value=t-f;
            bs=(boot(:,2,1,q)-boot(:,1,1,q))-(boot(:,2,2,q)-boot(:,1,2,q));
        end
        response=[response;interval_row("axis shift",metric(q),model,value, ...
            bs,nboot,seed)]; %#ok<AGROW>
    end
end
writetable(metrics,fullfile(here,'stripe_metrics.csv'));
writetable(intervals,fullfile(here,'stripe_intervals.csv'));
writetable(response,fullfile(here,'stripe_axis_response.csv'));
disp(metrics);
end

function m=stripe_metrics(xy,labels)
c=mean(xy,1); r=hypot(xy(:,1)-c(1),xy(:,2)-c(2));
centroids=zeros(5,2); width=zeros(5,1); fw90=width; angle=width;
for s=1:5
    points=xy(labels==s,:); assert(size(points,1)>20);
    centroids(s,:)=mean(points,1); p=points-centroids(s,:);
    sx=mean(p(:,1).^2); sy=mean(p(:,2).^2); covxy=mean(p(:,1).*p(:,2));
    theta=.5*atan2(2*covxy,sy-sx);
    normal=p(:,1)*cos(theta)-p(:,2)*sin(theta);
    width(s)=2.355*sqrt(mean(normal.^2));
    q=sample_quantile(normal,[.05 .95]); fw90(s)=q(2)-q(1);
    angle(s)=theta;
end
mean_angle=atan2(mean(sin(angle)),mean(cos(angle)));
projection=centroids*[cos(mean_angle);-sin(mean_angle)];
m=[c,sqrt(mean(r.^2)),2*sample_quantile(r,.99), ...
    mean(diff(sort(projection))),mean(width),mean(fw90),rad2deg(mean_angle)];
end

function q=sample_quantile(x,p)
x=sort(x(:)); pos=1+(numel(x)-1)*p;
q=x(floor(pos)).'+(pos-floor(pos)).*(x(ceil(pos)).'-x(floor(pos)).');
end

function r=interval_row(c,metric,model,value,bs,nboot,seed)
ci=sample_quantile(bs,[.025 .975]);
r=table(string(c),string(metric),string(model),value,ci(1),ci(2),nboot,seed, ...
    'VariableNames',{'case_name','metric','model','estimate','ci95_low', ...
    'ci95_high','n_bootstrap','seed'});
end
