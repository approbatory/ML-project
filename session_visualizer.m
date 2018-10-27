f = figure;
subplot(2,1,1);
p1 = plot(o.data.y.raw.full);

subplot(2,1,2);
[XL, ~, act,~,~,~,~,stats] = plsregress(o.data.X.fast, o.data.y.scaled, 2);
p2 = scatter(act(:,1), act(:,2), 4, o.data.y.scaled);


trial_start = find(diff(o.data.mask.fast) == 1);
trial_end = find(diff(o.data.mask.fast) == -1);

tr_ix = 100;
refresh(tr_ix, act, stats, o, trial_start, trial_end);


function search_for_point(ip, act, stats, o, trial_start, trial_end)
d = sum((ip(1:2) - act).^2,2);
[~, ind] = min(d);
qw = cumsum(o.data.mask.fast);
f_ind = find(qw == ind, 1);

tr_ids = cumsum(diff(o.data.mask.fast) == 1);
tr_ix = tr_ids(f_ind);

refresh(tr_ix, act, stats, o, trial_start, trial_end);
end

function trial_finder(ip, act, stats, o, trial_start, trial_end)
f_ind = round(ip(1));

tr_ids = cumsum(diff(o.data.mask.fast) == 1);
tr_ix = tr_ids(f_ind);
refresh(tr_ix, act, stats, o, trial_start, trial_end);
end

function refresh(tr_ix, act, stats, o, trial_start, trial_end)
tr_range = trial_start(tr_ix):trial_end(tr_ix);
subplot(2,1,1);
xl_ = xlim;
pp = plot(o.data.y.raw.full);
pp.ButtonDownFcn = @(src, event) trial_finder(event.IntersectionPoint, act, stats, o, trial_start, trial_end);
hold on;
plot(tr_range, o.data.y.raw.full(tr_range)); hold off;
xlim(xl_);
subplot(2,1,2);
sp = scatter(act(:,1), act(:,2), 20, o.data.y.raw.fast,...
    'filled', 'MarkerFaceAlpha', 0.05);
sp.ButtonDownFcn = @(src, event) search_for_point(event.IntersectionPoint, act, stats, o, trial_start, trial_end);
hold on;
xmean = mean(o.data.X.fast);
dat = (o.data.X.full(tr_range,:) - xmean) * stats.W;
scatter(dat(:,1), dat(:,2), 50, o.data.y.raw.full(tr_range), 'filled');
hold off;
end