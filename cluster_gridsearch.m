% add includes
clear all
addpath(genpath(pwd));
initdir = pwd;
basedir = '/projects/colin/SPINS_hcp/';
cd(basedir)

% options
link_fxn = 'ward';
dist_fxns = {'euclidean', 'seuclidean'};
clusters = [2:20];
iterations = 1000;

% load atlas
%cr = load_nii('/projects/colin/scratch2/epitome/hcp2/nii2/crad200.hcp.nii');

s=dir('SPN*');
for pdx = 1:size(s,1);
    subsy{pdx} = s(pdx).name;
end

% load t statistics
n=1;
for pdx = 1:size(s,1)
    try
        cd([basedir subsy{pdx}  '/scaled'])
        d = load_nii('spmT_0004.nii');
        data(n,:) = reshape(d.img,32767*3,1);
        name(n) = subsy(pdx);
        n=n+1;
    catch
        disp(sprintf([basedir subsy{pdx}  '/scaled/spmT_0004.nii not found']))
        clear d
    end
end
name=name';
data = data(:,1:32767*2);

% reshape atlas
%crad_rois = reshape(cr.img,32767*3,1);
%crad_rois = crad_rois(1:32767*2);
%for pdx = 1:size(data,2)
%    for rdx = 1:200
%        crad_dat(pdx,rdx) = mean(data(crad_rois==rdx,pdx));
%    end
%end

%% determine the number of clusters in the data
% Stability-Based Validation of Clustering Solutions Lange et al 2004

% holds outputs of grid search
gridsearch = zeros(length(dist_fxns), length(clusters));

i = 1; % distance loop
for dist_fxn = dist_fxns;

    j = 1; % n clusters loop
    for clst = clusters;

        distances = NaN([iterations,1]); % used for mean instability measure

        for iter = 1:iterations;

            % randomly split subjects into 2 groups
            idx = randperm(size(data, 1));
            idxa = idx(1:floor(length(idx)/2));
            idxb = idx(floor(length(idx)/2)+1:end);

            da = data(idxa, :);
            db = data(idxb, :);

            % cluster 'a' half of data
            Ya = pdist(da, dist_fxn{1});
            Za = linkage(squareform(Ya), link_fxn);
            Ca = cluster(Za,'maxclust', clst);

            % train knn classifier on this half of the data
            model = fitcknn(da, Ca, 'NumNeighbors', 4);

            % evaluate model (used for testing)
            %rsubloss = resubLoss(model);
            %model_cv = crossval(model);
            %kloss = kfoldLoss(model_cv);
            %disp(sprintf('model quality: resubstitution loss/kfold = %f/%f', rsubloss, kloss))

            % cluster 'b' half of data
            Yb = pdist(db, dist_fxn{1});
            Zb = linkage(squareform(Yb), link_fxn);
            Cb = cluster(Zb,'maxclust', clst);

            % predict cluster labels on b based on those found on a
            Cp = predict(model, db);

            % match labels across predicited partition (Cp, from Ca) and Cb
            % note we introduce all possible labels into both Cp and Cb to
            % allow the algorithm to function when some of the labels aren't
            % predicted (this happens if there are extremely small clusters)
            % in the data... typically due to a bad linkage function
            mappings = match_labels([Cp; unique(Cb)], [Cb; unique(Cb)]);
            Cpout = zeros(length(Cp), 1);
            for x = 1:clst;
                idx = find(Cp == x);
                Cpout(idx) = mappings(x);
            end
            Cp = Cpout; clear Cpout;

            % calculate cluster similarity, normalize by the theoretical
            % cluster distance obtained by chance. using a pessimistic null
            % probability, but it is easy to compute for now

            % distance = part_agree_coef(Cp, Cb); % rand index, depricated
            distance = pdist([Cp'; Cb'], 'hamming');
            distances(iter) = distance / (1-(1/clst));
        end

        % take mean of all successes and insert into grid search
        gridsearch(i,j) = nanmean(distances);
        disp(sprintf('distfxn=%i/%i, clst=%i/%i instability=%f', ...
                      i, length(dist_fxns), j, length(clusters), gridsearch(i,j)));

        j = j + 1; % n clusters loop
    end
    i = i + 1; % linkage loop
end

% plot findings
figure;
plot(1:length(clusters), squeeze(gridsearch(1,:)), 'color', 'black', 'linewidth', 2);
hold on;
plot(1:length(clusters), squeeze(gridsearch(2,:)), 'color', 'red', 'linewidth', 2);
ylabel('normalized instability');
xlabel('number of clusters');
set(gca, 'XTick', 1:length(clusters));
set(gca, 'XTickLabel', clusters);
box off;
legend(dist_fxns{1}, dist_fxns{2})
savefig('cluster_stability.fig')
clf;

% find the maxima, or the optimal values
optimum = min(min(min(gridsearch)));
[x,y] = ind2sub(size(gridsearch),find(gridsearch == optimum));

% save environment
save([initdir '/cluster_grid_search.mat'])

