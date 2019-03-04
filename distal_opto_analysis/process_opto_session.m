function [n_inhib, n_disinhib, n_total] = process_opto_session(load_path)
p_ = @(p) fullfile(load_path,p);

if ~exist(p_('opto.mat'), 'file')
    n_inhib = -2; n_disinhib = -2; n_total = -2; return;
end
load(p_('opto.mat'));

if exist(p_('ls'), 'dir') ~= 0
    ds = DaySummary(p_('distalopto.txt'), p_('ls'));
elseif exist(p_('clean'), 'dir') ~= 0
    ds = DaySummary(p_('distalopto.txt'), p_('clean'));
elseif exist(p_('proj'), 'dir') ~= 0
    ds = DaySummary(p_('distalopto.txt'), p_('proj'));
else
    error('no cells directory found');
end

sham_session = isfield(trial_inds, 'sham') && ~isempty(trial_inds.sham);
if sham_session
    fprintf('This is a sham session\n');
    sham_trials = trial_inds.sham;
    real_trials = trial_inds.real;
    off_trials = trial_inds.off;
else
    off_trials = trial_inds.off;
end

num_cells = ds.num_classified_cells;
num_trials = ds.num_trials;

%select traces during CS start to the end
%segment_length = 180;
cs_start = ds.trial_indices(:,2) - ds.trial_indices(:,1)+1;
%tr_end = cs_start + segment_length - 1; %only using 180 frames
tr_end = ds.trial_indices(:,4) - ds.trial_indices(:,1)+1;

res = detect_all_events(ds);
X_full_bin = Utils.binarize_res(res, ds.trial_indices(end));


%X = zeros(num_cells, segment_length, num_trials);
%X_b = X;
X = cell(1, num_trials); X_mean = zeros(num_cells, num_trials);
X_b = cell(1, num_trials); X_b_mean = zeros(num_cells, num_trials);
for tr_ix = 1:num_trials
    trial_length = size(ds.trials(tr_ix).traces,2);
    if (cs_start(tr_ix) > trial_length) || (tr_end(tr_ix) > trial_length)
        n_inhib = -1; n_disinhib = -1; n_total = -1;
        warning('%s: trial length too small: trial %d, length %d', load_path, tr_ix, trial_length);
        return;
    end
    X{tr_ix}(:,:) = ds.trials(tr_ix).traces(:, cs_start(tr_ix):tr_end(tr_ix));
    X_mean(:, tr_ix) = mean(X{tr_ix},2);
    X_b{tr_ix}(:,:) = X_full_bin(ds.trial_indices(tr_ix,2):ds.trial_indices(tr_ix,4),:).';
    X_b_mean(:,tr_ix) = mean(X_b{tr_ix},2);
end
%X_b = X_b ~= 0;

%mean_act = squeeze(sum(X_b, 2));
%mean_act = cellfun(@(x) mean(x,2) , X_b);
mean_act = X_b_mean;

if sham_session
    laser_sham = false(1, num_trials);
    laser_sham(sham_trials) = true;
    laser_real = false(1, num_trials);
    laser_real(real_trials) = true;
    laser_off = false(1, num_trials);
    laser_off(off_trials) = true;
    %mean_act_on = mean_act(:, laser_real);
    %mean_act_off = mean_act(:, laser_sham);
    mean_act_on = mean_act(:, laser_sham);
    mean_act_off = mean_act(:, laser_off);
else
    laser_on = true(1,num_trials);
    laser_on(off_trials) = false;
    mean_act_on = mean_act(:, laser_on);
    mean_act_off = mean_act(:, ~laser_on);
end



p = -ones(1,num_cells);
med_on = zeros(1,num_cells);
med_off = zeros(1,num_cells);
for ix = 1:num_cells
    p(ix) = ranksum(mean_act_on(ix, :), mean_act_off(ix, :));
    med_on(ix) = mean(mean_act_on(ix, :));
    med_off(ix)= mean(mean_act_off(ix,:));
end
fdr_BH = mafdr(p, 'BHFDR', true);
signif = fdr_BH < 0.1;

%median_mean_activity = median(mean_act.');
%sort by med_on - med_off to see split
%[~, off_ord] = sort(median_mean_activity);
%[~, off_ord] = sort(med_on - med_off);
%figure;
%plot(med_off(off_ord), '.k');
%hold on;
%scatter(1:num_cells, med_on(off_ord), 4, signif(off_ord));


%figure;
%boxplot(mean_act_off(off_ord,:).', 'plotstyle', 'compact');
%hold on;
%boxplot(mean_act_on(off_ord,:).', 'plotstyle', 'compact', 'colors', 'r');
%signif_indices = find(signif(off_ord));
%text(signif_indices, 6 + 0.*signif_indices, '*');


n_disinhib = sum(signif & (med_on > med_off));
n_inhib = sum(signif & (med_on < med_off));
n_total = sum(signif);

end