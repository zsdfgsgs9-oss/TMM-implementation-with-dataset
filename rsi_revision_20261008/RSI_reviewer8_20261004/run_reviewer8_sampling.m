function run_reviewer8_sampling()
% Sample the archived TMM initializer; retain bin counts, not 50M particles.
here = fileparts(mfilename('fullpath')); root = fileparts(here);
code = fullfile(root,'RSI大修仿真统一归档_20260911', ...
    '05_可复现脚本','root_matlab');
addpath(code,'-begin');
assert(strcmp(which('init_particles'),fullfile(code,'init_particles.m')));
cfg = struct('pattern_mode','spot','spot_radius',18e-9, ...
    'spot_center_x',0,'spot_center_y',0,'U0_max',.5,'alpha_max_deg',10.7);
seed=42; chunk_size=1e6; ns=[50000 50000000];
edges={linspace(0,.5,51),linspace(0,10.7,51)};
fine_edges=linspace(0,1,100001);
bins=table(); gof=table();
variables=["U0","alpha"]; units=["eV","deg"];
for count=ns
    rng(seed,'twister'); timer=tic; done=0;
    counts=zeros(2,50); fine_counts=zeros(2,100000);
    sums=zeros(2,1); sums2=sums;
    while done<count
        n=min(chunk_size,count-done);
        [v,p]=init_particles(n,cfg); %#ok<ASGLU>
        values={p.U0,rad2deg(p.alpha)};
        uniforms={p.U0/.5,sin(p.alpha).^2/sin(deg2rad(10.7))^2};
        for j=1:2
            counts(j,:)=counts(j,:)+histcounts(values{j},edges{j});
            fine_counts(j,:)=fine_counts(j,:)+histcounts(uniforms{j},fine_edges);
            sums(j)=sums(j)+sum(values{j});
            sums2(j)=sums2(j)+sum(values{j}.^2);
        end
        if count==50000
            exact_d=zeros(2,1);
            for j=1:2
                u=sort(uniforms{j}(:));
                exact_d(j)=max([(1:n)'/n-u;u-(0:n-1)'/n]);
            end
        end
        done=done+n;
        clear v p values uniforms
        if mod(done,1e7)==0 || done==count
            fprintf('Sampling: %d/%d particles counted\n',done,count);
        end
    end
    runtime=toc(timer);
    for j=1:2
        assert(sum(counts(j,:))==count && sum(fine_counts(j,:))==count);
        w=diff(edges{j}); density=counts(j,:)./(count*w);
        f=theory_cdf(edges{j},j);
        bin_theory=diff(f)./w;
        probability=counts(j,:)/count; z=1.959963984540054;
        center=(probability+z^2/(2*count))/(1+z^2/count);
        half=z*sqrt(probability.*(1-probability)/count+z^2/(4*count^2)) ...
            /(1+z^2/count);
        lower=max(0,center-half)./w; upper=min(1,center+half)./w;
        bins=[bins;table(repmat(count,50,1),repmat(variables(j),50,1), ...
            edges{j}(1:end-1)',edges{j}(2:end)', ...
            counts(j,:)',density',bin_theory',lower',upper', ...
            repmat(units(j),50,1),'VariableNames', ...
            {'N','variable','bin_low','bin_high','count','empirical_pdf', ...
            'theory_bin_pdf','ci95_pdf_low','ci95_pdf_high','unit'})]; %#ok<AGROW>
        if count==50000
            d=exact_d(j); bound=0; method="exact_sorted";
        else
            cumulative=cumsum(fine_counts(j,:))/count;
            before=[0 cumulative(1:end-1)];
            d=max([abs(cumulative-fine_edges(2:end)), ...
                abs(before-fine_edges(1:end-1))]);
            bound=1e-5; method="streaming_CDF_grid";
        end
        meanv=sums(j)/count;
        stdv=sqrt((sums2(j)-sums(j)^2/count)/(count-1));
        gof=[gof;table(count,variables(j),seed,chunk_size,method,d,bound, ...
            sqrt(mean((density-bin_theory).^2)),meanv,stdv,runtime, ...
            'VariableNames',{'N','variable','seed','chunk_size','gof_method', ...
            'CDF_max_deviation','CDF_grid_bound','bin_average_pdf_RMSE', ...
            'sample_mean','sample_std','runtime_s'})]; %#ok<AGROW>
    end
end
old=readtable(fullfile(root,'RSI第四版结论数据报告归档_20260916', ...
    '04_分项数据','D2_D3_径向与GOF','input_distribution_gof.csv'), ...
    'TextType','string');
for k=1:height(gof)
    dataset="generated_50k"; if gof.N(k)==50000000, dataset="generated_50m"; end
    row=old(old.dataset==dataset & old.variable==gof.variable(k),:);
    assert(height(row)==1 && abs(gof.CDF_max_deviation(k)-row.KS_D)<1e-10, ...
        'Initializer or chunk ordering no longer matches the audited sampler.');
end
writetable(bins,fullfile(here,'sampling_histogram_bins.csv'));
writetable(gof,fullfile(here,'sampling_gof.csv'));
disp(gof(:,{'N','variable','CDF_max_deviation','CDF_grid_bound'}));
end

function f=theory_cdf(x,j)
if j==1
    f=x/.5;
else
    f=sin(deg2rad(x)).^2/sin(deg2rad(10.7))^2;
end
end
