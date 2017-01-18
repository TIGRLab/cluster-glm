% load grid search from previous run
clear all
addpath(genpath(pwd));
basedir = '/projects/colin/SPINS_hcp/';
load('cluster_grid_search.mat');
load('includes/SPM_bat_scripts/group_t12.mat'); % load base SPM struct for running GLMs

%% options
% CON_NAMES is in the same order as the output from SPM, DO NOT REORDER
CON_NAMES={'OBS_EM_NE'; 'OBS_EM_FX'; 'IMI_EM_NE'; 'IMI_EM_FX'; 'EM_IMI_OBS'; 'NE_IMI_OBS'; ...
    'OBS_FE_FX';  'OBS_SA_FX';  'OBS_HA_FX';  'OBS_AN_FX'; ...
    'IMI_FE_FX';  'IMI_SA_FX';  'IMI_HA_FX';  'IMI_AN_FX'; ...
    'OBS_FE_NE';  'OBS_SA_NE';  'OBS_HA_NE';  'OBS_AN_NE'; ...
    'IMI_FE_NE';  'IMI_SA_NE';  'IMI_HA_NE';  'IMI_AN_NE'; ...
    'main_OBS_NE'; 'main_OBS_FE'; 'main_OBS_SA'; 'main_OBS_HA'; 'main_OBS_AN'; 'main_OBS_FX'; ...
    'main_IMI_NE'; 'main_IMI_FE'; 'main_IMI_SA'; 'main_IMI_HA'; 'main_IMI_AN'; 'main_IMI_FX'; ...
};

con = 5; % index from above list
ana = 'scaled'; % input filename prefix
clusname = 'all_3solution_euclidean_ward'; % output directory name

%% plot standard deviation of t-stats across each subject
figure; hist(std(data')); title('t-statistics'); savefig('stat_hist.fig'); clf;

%% plot the best clustering solution found
dist_fxn = dist_fxns{x};
link_fxn = 'ward';
n_clust = clusters(z);

Y = pdist(data, dist_fxn);
Z = linkage(squareform(Y), link_fxn);
C = cluster(Z,'maxclust', n_clust);

% find the height that would dendrogram into n clusters
thresh = mean(Z(end-(n_clust-2:n_clust-1), 3));

figure;
subplot(1,2,1);
[h t idx] = dendrogram(Z, 0, 'Orientation','left', 'ColorThreshold', thresh);
subplot(1,2,2);

display = squareform(Y);
%imagesc(flipud(data(idx, :)), [-3 3]);
imagesc(flipud(display(idx, idx)));
colormap(hot);
colorbar;

%% export 2nd-level GLM of each cluster
con_name = CON_NAMES{con};

for clus = 1:max(C)
    condir = [basedir clusname num2str(clus) '/' ];
    mkdir(condir)

    job.dir = {[condir '/' con_name  '/']};
    mkdir([condir '/' con_name  '/'])
    cd([ condir '/' con_name  '/']);

    %clear scans in current batch
    job.des.t1.scans = {};

    n = 1;
    for pdx = 1:length(C)
        if C(pdx) == clus
            job.des.t1.scans(n) = {[basedir name{pdx}  '/' ana '/' sprintf('con_%04d', con) '.nii,1']};
            n = n + 1;
        end
    end

    try
        if ~isempty(ls('SPM.mat'))
            delete SPM.mat
        end;
    end

    % spm settings
    job.masking.tm.tm_none=1;
    job.masking.im= 0;
    job.masking.em= {''};
    job.globalc.g_omit=1;
    job.globalm.glonorm=1;
    job.globalm.gmsca.gmsca_no=1;

    spm_run_factorial_design(job) % generate SPM design matrix
    load SPM.mat                  % load SPM design matrix
    SPM=spm_spm(SPM)              % estimate design matrix
    save('SPM', 'SPM')

    cons = [1; -1];
    names = {'pos', 'neg'};
    curdir = pwd;
    analyze_spm_contrasts(curdir, cons, names);

end

save([initdir '/cluster_final.mat'])

