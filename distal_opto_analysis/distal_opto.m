%cd _distalopto_traces/

%cd enphr/

%cd m953/

%cd m953-0122/


load_path = fullfile('~','_distalopto_traces', 'mcherry',...
    'm892', 'm892-1016');
p_ = @(p) fullfile(load_path,p);

load(p_('opto.mat'));


ds = DaySummary(p_('distalopto.txt'), p_('ls'));

off_trials = trial_inds.off;

num_cells = ds.num_classified_cells;
num_trials = ds.num_trials;
%%
%select traces during CS start to the end
segment_length = 180;
cs_start = ds.trial_indices(:,2) - ds.trial_indices(:,1)+1;
tr_end = cs_start + segment_length - 1; %only using 180 frames

res = detect_all_events(ds);
X_full_bin = Utils.binarize_res(res, ds.trial_indices(end));
%X_trials_bin = cell(1,num_trials);
%for tr_ix = 1:num_trials
%    X_trials_bin{tr_ix} = X_full_bin(ds.trial_indices(tr_ix,1):ds.trial_indices(tr_ix,4),:);
%end
%%
X = zeros(num_cells, segment_length, num_trials);
X_b = X;
for tr_ix = 1:num_trials
    X(:,:,tr_ix) = ds.trials(tr_ix).traces(:, cs_start(tr_ix):tr_end(tr_ix));
    X_b(:,:,tr_ix) = X_full_bin(ds.trial_indices(tr_ix,2):ds.trial_indices(tr_ix,2)+segment_length-1,:).';
end
X_b = X_b ~= 0;
%%
mean_act = squeeze(sum(X_b, 2));

laser_on = true(1,num_trials);
laser_on(off_trials) = false;

mean_act_on = mean_act(:, laser_on);
mean_act_off = mean_act(:, ~laser_on);

p = -ones(1,num_cells);
med_on = zeros(1,num_cells);
med_off = zeros(1,num_cells);
for ix = 1:num_cells
    p(ix) = ranksum(mean_act_on(ix, :), mean_act_off(ix, :));
    %med_on(ix) = median(mean_act_on(ix, :));
    %med_off(ix)= median(mean_act_off(ix,:));
    med_on(ix) = mean(mean_act_on(ix, :));
    med_off(ix)= mean(mean_act_off(ix,:));
end
fdr_BH = mafdr(p, 'BHFDR', true);
signif = fdr_BH < 0.1;

median_mean_activity = median(mean_act.');
%sort by med_on - med_off to see split
%[~, off_ord] = sort(median_mean_activity);
[~, off_ord] = sort(med_on - med_off);
%figure;
%plot(med_off(off_ord), '.k');
%hold on;
%scatter(1:num_cells, med_on(off_ord), 4, signif(off_ord));


figure;
boxplot(mean_act_off(off_ord,:).', 'plotstyle', 'compact');
hold on;
boxplot(mean_act_on(off_ord,:).', 'plotstyle', 'compact', 'colors', 'r');
signif_indices = find(signif(off_ord));
text(signif_indices, 6 + 0.*signif_indices, '*');

%%
num_excited = sum(signif & (med_on > med_off));
num_inhibited = sum(signif & (med_on < med_off));