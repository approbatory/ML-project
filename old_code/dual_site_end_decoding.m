ds_hpc = quick_ds('../cohort14_dual/c14m6/c14m6d10', 'deprobe', 'nocells', 'cm', 'hpc_cm01_fix');
ds_prl = quick_ds('../cohort14_dual/c14m6/c14m6d10', 'deprobe', 'nocells', 'cm', 'prl_cm01_fix');
dayset = my_daysets('c14m6');
dayset = dayset(2);

points = 0:0.1:1;
algs = my_algs({'mvnb2', 'linsvm'}, {'original'});

[tr_hpc, te_hpc] = end_dec(ds_hpc, algs, points);
[tr_prl, te_prl] = end_dec(ds_prl, algs, points);
[tr_both, te_both] = end_dec([ds_hpc, ds_prl], algs, points);

figure;
subplot(1,3,1);
plotter(tr_hpc, te_hpc, algs, points, 'hpc');
subplot(1,3,2);
plotter(tr_prl, te_prl, algs, points, 'prl');
subplot(1,3,3);
plotter(tr_both, te_both, algs, points, 'both');

%%
algs = my_algs({'mvnb2', 'ecoclin'}, {'original'});
[tr_place(1,:), te_place(1,:), s_tr_place(1,:), s_te_place(1,:), m_tr_place(1,:), m_te_place(1,:)] = place_dec(ds_hpc, algs);
[tr_place(2,:), te_place(2,:), s_tr_place(2,:), s_te_place(2,:), m_tr_place(2,:), m_te_place(2,:)] = place_dec(ds_prl, algs);
[tr_place(3,:), te_place(3,:), s_tr_place(3,:), s_te_place(3,:), m_tr_place(3,:), m_te_place(3,:)] = place_dec([ds_hpc, ds_prl], algs);

d(1).label = 'hpc';
d(2).label = 'prl';
d(3).label = 'both';
%% comparing filling
fils(1).label = 'binary';
fils(2).label = 'box';
fils(3).label = 'copy';

for i = 1:numel(fils)
    [tr_place_fil(i,:), te_place_fil(i,:), s_tr_place_fil(i,:), s_te_place_fil(i,:), m_tr_place_fil(i,:), m_te_place_fil(i,:)] = place_dec([ds_hpc, ds_prl], algs, fils(i).label);
end
%%
figure;
subplot(3,1,1);
plotmat(sprintf('Place decoding: %s, both', dayset.label), tr_place_fil, te_place_fil, algs, fils, 1);
subplot(3,1,2);
plotmat(sprintf('Place decoding: %s, both, tested on moving', dayset.label), s_tr_place_fil, s_te_place_fil, algs, fils, 1);
subplot(3,1,3);
plotmat(sprintf('Place decoding: %s, both, only moving', dayset.label), m_tr_place_fil, m_te_place_fil, algs, fils, 1);
%%
figure;
subplot(3,1,1);
plotmat(sprintf('Place decoding: %s', dayset.label), tr_place, te_place, algs, d, 2);
subplot(3,1,2);
plotmat(sprintf('Place decoding: %s, tested on moving', dayset.label), s_tr_place, s_te_place, algs, d, 2);
subplot(3,1,3);
plotmat(sprintf('Place decoding: %s, only moving', dayset.label), m_tr_place, m_te_place, algs, d, 2);
%%
function [tr, te, s_tr, s_te, m_tr, m_te] = place_dec(ds, algs, fil)
if ~exist('fil', 'var')
    fil = 'binary';
end
PARLOOPS = 64;

subset = subset_moving(ds(1));

X = cell(1,numel(ds)); ks = cell(1,numel(ds));
for ix = 1:numel(ds)
    [X{1,ix}, ks{1,ix}, errf] = ds_dataset(ds(ix), 'filling', fil);
    if ~isequal(ks{1,ix}, ks{1,1})
        error('unequal ks');
    end
end
X = cell2mat(X);
ks = ks{1,1};
for i = 1:numel(algs)
    [tr{i}, te{i}, s_tr{i}, s_te{i}, m_tr{i}, m_te{i}] = evaluate_alg(algs(i), X, ks, 'eval_f', errf,...
        'train_frac', 0.7, 'par_loops', PARLOOPS, 'subset', subset);
end
end
%%

function [tr, te] = end_dec(ds, algs, points)
PARLOOPS = 512;

trials_chosen = strcmp({ds(1).trials.start}, 'west'); % using the changing path (west for both)
indices_chosen = 1:150;
indices_chosen = indices_chosen(trials_chosen);


tr = cell(numel(algs), numel(points)); te = tr; tr_sh = tr; te_sh = tr; 
for j = 1:length(points)
    X_ = cell(1,numel(ds)); ks_ = cell(1,numel(ds));
    for ix = 1:numel(ds)
        [X_{1,ix}, ks_{1,ix}, em_] = ds_dataset(ds(ix), 'selection', points(j), 'filling', 'binary',...
            'trials', trials_chosen, 'target', {ds(ix).trials.end});
        if ~isequal(ks_{1,ix}, ks_{1,1})
            error('unequal ks_');
        end
    end
    X_ = cell2mat(X_);
    ks_ = ks_{1,1};
    for i = 1:numel(algs)
        alg = algs(i);
        [tr{i,j}, te{i,j}] = evaluate_alg(alg, X_, strcmp(ks_,'north'),...
            'eval_f', @(k,p) mean(~(k(:)==p(:))), 'par_loops', PARLOOPS);
        %[tr_sh{i,j}, te_sh{i,j}] = evaluate_alg(alg, @() shuffle(X_,strcmp(ks_,'north')), strcmp(ks_,'north'),...
        %    'eval_f', @(k,p) mean(~(k(:)==p(:))), 'X_is_func', true, 'par_loops', PARLOOPS);
    end
end
end


function plotter(tr, te, algs, points, region)

tr_mean = cellfun(@mean, tr);
te_mean = cellfun(@mean, te);

%tr_sh_mean = cellfun(@mean, tr_sh);
%te_sh_mean = cellfun(@mean, te_sh);


tr_un = cellfun(@(x) std(x)/sqrt(length(x)), tr);
te_un = cellfun(@(x) std(x)/sqrt(length(x)), te);

%tr_sh_un = cellfun(@(x) std(x)/sqrt(length(x)), tr_sh);
%te_sh_un = cellfun(@(x) std(x)/sqrt(length(x)), te_sh);


%colors = 'rgbcymk';
colors = 'rb';

%figure;
%subplot(1,2,1);
for i = 1:numel(algs)
    errorbar(points, tr_mean(i,:), tr_un(i,:), ['-.' colors(i)], 'DisplayName', [algs(i).name ' train']);
    hold on;
    errorbar(points, te_mean(i,:), te_un(i,:), ['-' colors(i)], 'DisplayName', [algs(i).name ' test']);
    hold on;
end
ylim([0 0.3]);
legend location best
xlabel('Fraction of turn completed');
ylabel('End arm classification error');
title(sprintf('End arm prediction at positions along turn (%s)', region));

%subplot(1,2,2);
%for i = 1:numel(algs)
%    errorbar(points, tr_sh_mean(i,:), tr_sh_un(i,:), ['-.' colors(i)], 'DisplayName', [algs(i).name ' train']);
%    hold on;
%    errorbar(points, te_sh_mean(i,:), te_sh_un(i,:), ['-' colors(i)], 'DisplayName', [algs(i).name ' test']);
%    hold on;
%end
%ylim([0 0.3]);
%legend location best
%xlabel('Fraction of turn completed');
%ylabel('End arm classification error');
%title(sprintf('End arm prediction at positions along turn (shuffles of %s)', region));
end