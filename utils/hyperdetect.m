function events = hyperdetect(traces, varargin)
p = inputParser;
p.addRequired('traces', @(x) size(x,1)>size(x,2));
p.addParameter('z_val', 3, @isscalar);
p.addParameter('Progress', false, @islogical);
p.addParameter('OutType', 'full', @(x) ~isempty(validatestring(x, {'full', 'onset', 'peak', 'mid', 'predip'})));
p.parse(traces, varargin{:});
opt = p.Results;

[~, n_cells] = size(traces);
events = zeros(size(traces));

if opt.Progress
    %progressbar('cells...');
    fprintf('Detecting %d cells\n', n_cells);
end

for c_ix = 1:n_cells
    events(:,c_ix) = hyperdetect_onecell(traces(:,c_ix), opt.OutType, opt.z_val);
    if opt.Progress
        %progressbar(c_ix/n_cells);
        if mod(c_ix, 200) == 0
            fprintf('%d\n', c_ix);
        end
    end
end
if opt.Progress
    fprintf('Done\n');
end
end


function events = hyperdetect_onecell(tr, out_type, z_val)
events = zeros(size(tr));

tr_dn = wdenoise(double(tr));
noise = tr - tr_dn;

sigma = std(noise);

min_mag = z_val * sigma;

peak_th = prctile(tr_dn, 25) + min_mag;

[vals, idx] = findpeaks(tr_dn);

%filter by peak_th
high_peaks = vals > peak_th;

vals = vals(high_peaks);
idx = idx(high_peaks);

is_decreasing = [false; diff(tr_dn) <= 0];

[base_vals, peak_start, th_start] = deal(zeros(size(vals)));

for p_i = 1:numel(idx)
    peak_start(p_i) = idx(p_i) - 1;
    while peak_start(p_i) > 1 && ~is_decreasing(peak_start(p_i))
        peak_start(p_i) = peak_start(p_i) - 1;
    end
    if peak_start(p_i) < 1
        peak_start(p_i) = 1;
    end
    base_vals(p_i) = tr_dn(peak_start(p_i));
    
    th_start(p_i) = idx(p_i) - 1;
    while th_start(p_i) > 1 &&...
            tr_dn(th_start(p_i)) - base_vals(p_i) > min_mag &&...
            ~is_decreasing(th_start(p_i))
        th_start(p_i) = th_start(p_i) - 1;
    end
    if th_start(p_i) < 1
        th_start(p_i) = 1;
    end
end
increase_vals = vals - base_vals;

keep = increase_vals > min_mag;

increase_vals = increase_vals(keep);
peak_start = peak_start(keep);
th_start = th_start(keep);
idx = idx(keep);

for j = 1:numel(th_start)
    switch out_type
        case 'full'
            events(th_start(j):idx(j)) = increase_vals(j);
        case 'onset'
            events(th_start(j)) = increase_vals(j);
        case 'peak'
            events(idx(j)) = increase_vals(j);
        case 'mid'
            midpoint = round((th_start(j) + idx(j))/2);
            events(midpoint) = increase_vals(j);
        case 'predip'
            events(peak_start(j):idx(j)) = increase_vals(j);
    end
end
end