%% script parameters
visualize = false; %use to make intermediate plots

%% loading file

data_options = { 'Mouse-2024-20150315_090450-linear-track-TracesAndEvents.mat',...
    'Mouse-2024-20150317_202708-linear-track-TracesAndEvents.mat'};
chosen_option = data_options{2};

path_to_data = fullfile('../linear_track', chosen_option);

data_struct = load(path_to_data);
tracesEvents = data_struct.tracesEvents;

%% constructing dataset
full_X = tracesEvents.rawTraces;
full_y = tracesEvents.position;

[forward_ranges, backward_ranges] = select_directions(full_y, 'ranges', true);

%% TODO
n_forward_trials = numel(forward_ranges.start);
trials_X = cell(n_forward_trials,1);
trials_y = cell(n_forward_trials,1);
for i = 1:n_forward_trials
    s = forward_ranges.start(i);
    e = forward_ranges.end(i);
    trials_X{i} = full_X(s:e,:);
    trials_y{i} = full_y(s:e,1);
end

%% plot verification
if visualize
    figure; hold on;
    n_count = 0;
    for i = 1:n_forward_trials
        plot(n_count + (1:length(trials_y{i})),trials_y{i}(:,1), '-o');
        n_count = n_count + length(trials_y{i});
    end
end
%% PLS on 1/3 of data:
%[pls_X, pls_y, pls_sel] = subsel(1/3, trials_X, trials_y);
pls_X = cell2mat(trials_X);
pls_y = cell2mat(trials_y);
pls_y = 100.*(pls_y - min(pls_y))./range(pls_y); %100cm range
max_dim = 100;
[L, ~, ~, ~, ~, ~, MSE] = plsregress(pls_X, pls_y, max_dim, 'cv', 10);
[L_n, ~, ~, ~, ~, ~, MSE_n] = plsregress(pls_X, pls_y, max_dim, 'cv', 'resubstitution');

figure; hold on;
plot(0:max_dim, sqrt(MSE(2,:)), '-o');
plot(0:max_dim, sqrt(MSE_n(2,:)), '-o');
legend 'test error' 'training error'
xlabel('PLS dimensions');
ylabel('RMS error (cm)');
title('RMS position error per PLS dimension');
%% my own CV
cat_y = cell2mat(trials_y);
rescaler = @(y) 100.*(y - min(cat_y))./range(cat_y);

max_dim = 20;
max_reps = 400;
RMSE_n = zeros(max_reps, max_dim);
RMSE = zeros(max_reps, max_dim);
for dim = 1:max_dim
    parfor rep = 1:max_reps
        [pX, py, p_sel] = subsel(0.7, trials_X, trials_y); pX = cell2mat(pX); py = cell2mat(py);
        py = rescaler(py);
        
        [XL, YL, XS, YS, BETA] = plsregress(pX, py, dim);
        predy_n = [ones(length(py),1), pX]*BETA;
        RMSE_n(rep,dim) = sqrt(mean((py - predy_n).^2));
        [testX, testy] = subsel(~p_sel, trials_X, trials_y); testX = cell2mat(testX); testy = cell2mat(testy);
        testy = rescaler(testy);
        predy = [ones(length(testy),1), testX]*BETA;
        RMSE(rep,dim) = sqrt(mean((testy - predy).^2));
    end
    fprintf('%d ', dim); if mod(dim,10)==0, fprintf('\n'); end
end
fprintf('\n');
%%
m_RMSE = mean(RMSE);
e_RMSE = std(RMSE)./sqrt(max_reps);
m_RMSE_n = mean(RMSE_n);
e_RMSE_n = std(RMSE_n)./sqrt(max_reps);
figure; hold on;
errorbar(m_RMSE, e_RMSE);
errorbar(m_RMSE_n, e_RMSE_n);
legend test train
xlabel 'PLS dimensions'
ylabel 'RMS error (cm)'
title 'PLS train and test RMS error vs. reduced space dimension'

%% using 5 PLS dimensions
PLS_DIMS = 30;

[pX, py, p_sel] = subsel(0.5, trials_X, trials_y); pX = cell2mat(pX); py = cell2mat(py);
py = rescaler(py);
[dX, dy, d_sel] = subsel(~p_sel, trials_X, trials_y); dX = cell2mat(dX); dy = cell2mat(dy);
dy = rescaler(dy);

[XL, yL, XS] = plsregress(pX, py, PLS_DIMS);
%reconstructed as XS_recon = pX/(XL');
XS_test = dX/(XL');

figure;
subplot(2,5,1); scatter3(XS(:,1), XS(:,2), XS(:,3), 5, py);
subplot(2,5,2); scatter3(XS_test(:,1), XS_test(:,2), XS_test(:,3), 5, dy);

%TODO calculate (d')^2 in XS_test, for the various spatial bins
%produce the minimal d', as well as 1-to-others for each bin
%for number of neurons less than 5, use no PLS, otherwise use PLS dim 5
num_bins = 20;
midLen = 100;
my_binner = @(y) gen_place_bins(y, num_bins, midLen);
dks = my_binner(dy);
X_cell = cell(num_bins,1);
for i = 1:num_bins
    X_cell{i} = XS_test(dks==i,:);
end

dp2_matrix = dprime2_mat(X_cell, false);
dp2_other = dprime2_other(X_cell, false);

subplot(2,5,3); imagesc(dp2_matrix, [0 120]); colorbar; title 'pairwise (d'')^2'
subplot(2,5,4); plot(dp2_other, '-o'); title 'one to rest (d'')^2'
subplot(2,5,5); plot(min(dp2_matrix), '-o'); title 'nearest (d'')^2'
%
dX_shuf = shuffle(dX, dks);
dX_cell = cell(num_bins,1);
dX_cell_shuf = cell(num_bins,1);
XS_shuf_test = dX_shuf/(XL');
X_cell_shuf = cell(num_bins,1);
for i = 1:num_bins
    X_cell_shuf{i} = XS_shuf_test(dks==i,:);
    dX_cell{i} = dX(dks==i,:);
    dX_cell_shuf{i} = dX_shuf(dks==i,:);
end
dp2_matrix_nopls = dprime2_mat(dX_cell, true);
dp2_other_nopls = dprime2_other(dX_cell, true);

dp2_matrix_shuf = dprime2_mat(X_cell_shuf, false);
dp2_other_shuf = dprime2_other(X_cell_shuf, false);

dp2_matrix_shuf_nopls = dprime2_mat(dX_cell_shuf, true);
dp2_other_shuf_nopls = dprime2_other(dX_cell_shuf, true);

subplot(2,5,6); scatter3(XS(:,1), XS(:,2), XS(:,3), 5, py);
subplot(2,5,7); scatter3(XS_shuf_test(:,1), XS_shuf_test(:,2), XS_shuf_test(:,3), 5, dks);
subplot(2,5,8); imagesc(dp2_matrix_shuf, [0 120]); colorbar; title 'pairwise (d'')^2, on shuffle'
subplot(2,5,4); hold on; plot(dp2_other_shuf, '-o'); %plot(dp2_other_nopls, '-o'); plot(dp2_other_shuf_nopls, '-o'); 
legend unshuffled shuffled 'unshuffled diag' 'shuffled diag'
subplot(2,5,5); hold on; plot(min(dp2_matrix_shuf), '-o'); %plot(min(dp2_matrix_nopls), '-o'); plot(min(dp2_matrix_shuf_nopls), '-o'); 
legend unshuffled shuffled 'unshuffled diag' 'shuffled diag'
subplot(2,5,9); imagesc(dp2_matrix_nopls, [0 120]); colorbar; title 'pairwise (d'')^2, diag unshuffled'
subplot(2,5,10); imagesc(dp2_matrix_shuf_nopls, [0 120]); colorbar; title 'pairwise (d'')^2, diag shuffled'
%

%%
X_a = randn(1000, 10);
X_b = randn(1000, 10) + 1;
dp2 = dprime2_CV(X_a, X_b, false);
dp2_shuf = dprime2_CV(X_a, X_b, true);

%%

%% helper functions
function varargout = subsel(f, varargin)
L = cellfun(@length, varargin);
assert(~isempty(L), 'give at least one array');
assert(all(L == L(1)), 'must all be of equal length');
L = L(1);
if isscalar(f) && isnumeric(f)
    sel = randperm(L) <= L*f;
elseif islogical(f)
    sel = f;
end
varargout = cellfun(@(x) x(sel), varargin, 'UniformOutput', false);
varargout{length(varargin) + 1} = sel;
end

function [dp2_matrix_test, dp2_matrix_train] = dprime2_mat(X_cell, diag_only)
dp2_matrix_test = zeros(numel(X_cell));
dp2_matrix_train = zeros(numel(X_cell));
for i = 1:numel(X_cell)
    for j = 1:numel(X_cell)
        if i==j
            dp2_matrix_test(i,j) = Inf;
            dp2_matrix_train(i,j) = Inf;
        else
            [dp2_matrix_test(i,j), dp2_matrix_train(i,j)] = dprime2_CV(X_cell{i}, X_cell{j}, diag_only);
        end
    end
end
end

function [dp2_other_test, dp2_other_train] = dprime2_other(X_cell, diag_only)
dp2_other_test = zeros(1,numel(X_cell));
dp2_other_train = zeros(1,numel(X_cell));
for i = 1:numel(X_cell)
    X_a = X_cell{i};
    X_b = cell2mat(X_cell((1:numel(X_cell))~=i));
    [dp2_other_test(i), dp2_other_train(i)] = dprime2_CV(X_a, X_b, diag_only);
end
end

function dp2 = dprime2(X_a, X_b, diag_only)
mu_a = mean(X_a).';
mu_b = mean(X_b).';
delta_mu = mu_a - mu_b;
if ~diag_only
    Sigma2 = cov([(X_a - mu_a.'); (X_b - mu_b.')]);
    w_opt = Sigma2\delta_mu;
    dp2 = delta_mu.' * w_opt;
else
    Sigma2_diag = var([(X_a - mu_a.'); (X_b - mu_b.')]).';
    w_diag = Sigma2_diag.\delta_mu;
    dp2 = delta_mu.' * w_diag;
end
%w_diag = delta_mu./diag(Sigma2);

%dp_shuffle2 = delta_mu.' * w_diag;
%dp_diagonal2 = dp2^2 / (w_diag.' * Sigma2 * w_diag);
end


function [dp2_test, dp2_train] = dprime2_CV(X_a, X_b, diag_only, k_fold)
if ~exist('k_fold', 'var')
    k_fold = 2;
end
dp2_test = zeros(k_fold, 1);
dp2_train = zeros(k_fold, 1);
for k_i = 1:k_fold
    [X_a_train, X_a_test, X_b_train, X_b_test] = kfold_selector(k_fold, k_i, X_a, X_b);

    mu_a_train = mean(X_a_train).';
    mu_b_train = mean(X_b_train).';
    delta_mu_train = mu_a_train - mu_b_train;
    if ~diag_only
        Sigma2_train = cov([(X_a_train - mu_a_train.'); (X_b_train - mu_b_train.')]);
        w_opt_train = Sigma2_train\delta_mu_train;
    else
        Sigma2_diag_train = var([(X_a_train - mu_a_train.'); (X_b_train - mu_b_train.')]).';
        w_opt_train = Sigma2_diag_train.\delta_mu_train;
        dp2_train(k_i) = delta_mu_train.' * w_opt_train;
    end
    delta_mu_test = mean(X_a_test).' - mean(X_b_test).';
    dp2_test(k_i) = delta_mu_test.' * w_opt_train;
end
dp2_test = mean(dp2_test);
dp2_train = mean(dp2_train);
end