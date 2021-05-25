function events = wavelet_event_detection(traces, varargin)
p = inputParser;
p.addRequired('traces', @(x) size(x,1)>size(x,2));
p.addParameter('z_val', 4, @isscalar);
p.addParameter('fps', 25, @isscalar);
p.addParameter('MAX_IT', 10, @isscalar);
p.addParameter('med_window_duration', seconds(50), @isscalar);
p.addParameter('Progress', false, @islogical);
p.addParameter('OutType', 'full', @(x) ~isempty(validatestring(x, {'full', 'onset', 'peak', 'mid'})));
p.addParameter('KernelWidth', 1, @isscalar);
p.parse(traces, varargin{:});
opt = p.Results;

[n_frames, n_cells] = size(traces);
events = zeros(n_frames, n_cells);

%traces_dn = wdenoise(traces);
%med = median(traces_dn);

%traces_dn_med = traces_dn - med;

rstd = @(x) 1.4826*mad(x,1);

if opt.Progress
    progressbar('cells...');
end
for c_ix = 1:n_cells
    tr = traces(:,c_ix);
    tr_dn = wdenoise(tr);
    med = median(tr_dn);
    tr_dn_med = tr_dn - med;
    
    %tr_dn = traces_dn(:,c_ix);
    %tr_dn_med = traces_dn_med(:,c_ix);
    
    mu = median(tr_dn_med);
    sigma = rstd(tr_dn_med);
    
    peaks = tr_dn_med > mu + opt.z_val*sigma;
    shade = shadow(tr_dn_med, peaks);
    
    for i = 1:opt.MAX_IT
        tr_dn_nopeaks = tr_dn;
        tr_dn_nopeaks(shade) = nan;
        
        med_baseline = medfilt1(tr_dn_nopeaks, opt.fps * seconds(opt.med_window_duration), 'omitnan');
        med_baseline(shade) = nan;
        
        good_frames = find(~isnan(med_baseline));
        good_values = med_baseline(good_frames);
        med_baseline = interp1(good_frames, good_values, (1:numel(med_baseline)).');
        
        tr_dn_med = tr_dn - med_baseline;
        
        non_peaks = tr_dn_med(~shade);
        
        mu = median(non_peaks);
        sigma = rstd(non_peaks);
        
        peaks_new = tr_dn_med > mu + opt.z_val*sigma;
        shade_new = shadow(tr_dn_med, peaks_new);
        
        if isequal(peaks_new, peaks)
            peaks = peaks_new;
            shade = shade_new;
            exit_iters = i;
            break;
        end
        
        peaks = peaks_new;
        shade = shade_new;
    end
    
    peak_onsets = find(diff([0;peaks]) == 1);
    peak_offsets = find(diff([peaks;0]) == -1);
    
    assert(numel(peak_onsets) == numel(peak_offsets));
    assert(all(peak_offsets >= peak_onsets));
    
    is_decreasing = [false; diff(tr_dn) < 0];
    
    peak_magnitudes = [];
    peak_inits = [];
    peak_tops = [];
    
    for i = 1:numel(peak_onsets)
        on = peak_onsets(i);
        off = peak_offsets(i);
        
        segment = tr_dn(on:off); %inclusive
        
        if numel(segment) < 3
            skip_findpeaks = true;
        else
            [vals, idx] = findpeaks(segment);
            skip_findpeaks = isempty(vals);
        end
        if skip_findpeaks
            [vals, idx] = max(segment);
        end
        idx = idx + on - 1;
        
        [base_vals, peak_start] = deal(zeros(size(vals)));
        for p_i = 1:numel(idx)
            peak_start(p_i) = idx(p_i) - 1;
            while peak_start(p_i) > 1 && ~is_decreasing(peak_start(p_i))
                peak_start(p_i) = peak_start(p_i) - 1;
            end
            if peak_start(p_i) < 1
                peak_start(p_i) = 1;
            end
            base_vals(p_i) = tr_dn(peak_start(p_i));
        end
        increase_vals = vals - base_vals;
        
        peak_magnitudes = [peak_magnitudes ; increase_vals(:)];
        peak_inits = [peak_inits ; peak_start(:)];
        peak_tops = [peak_tops ; idx(:)];
    end
    
    %filter by mags > z*sigma
    keep = peak_magnitudes > opt.z_val*sigma;
    peak_magnitudes = peak_magnitudes(keep);
    peak_inits = peak_inits(keep);
    peak_tops = peak_tops(keep);
    
    %ev_trace = 0*tr_dn;
    for j = 1:numel(peak_inits)
        switch opt.OutType
            case 'full'
                events(peak_inits(j):peak_tops(j), c_ix) = peak_magnitudes(j);
            case 'onset'
                events(peak_inits(j), c_ix) = peak_magnitudes(j);
            case 'peak'
                events(peak_tops(j), c_ix) = peak_magnitudes(j);
            case 'mid'
                midpoint = round( (peak_inits(j) + peak_tops(j))/2 );
                events(midpoint, c_ix) = peak_magnitudes(j);
        end
    end
    
    %if opt.KernelWidth ~= 1
    %    ev_tmp = conv(events(:, c_ix), ones(opt.KernelWidth,1), 'full');
    %end
    
    if opt.Progress
        progressbar(c_ix/n_cells);
    end
end

if opt.KernelWidth > 1
    X_ = conv2(events, ones(opt.KernelWidth, 1), 'full');
    X_ = X_(1:size(events,1),:);
    events = X_;
end

end

function sh = shadow(t, pk)
t = t(:);
pk = pk(:);

assert(numel(t) == numel(pk), 'Trace and peak indicator do not match in size');

dt = diff(t);
is_decreasing = [false; dt < 0];
is_increasing = [false; dt > 0];

peak_onsets = find(diff([0;pk]) == 1);
peak_offsets = find(diff([pk;0]) == -1);

assert(numel(peak_onsets) == numel(peak_offsets));
assert(all(peak_offsets >= peak_onsets));

sh = false(size(pk));

for i = 1:numel(peak_onsets)
    p_on = peak_onsets(i);
    p_off = peak_offsets(i);
    
    sh_on = p_on - 1;
    while sh_on > 1 && is_increasing(sh_on)
        sh_on = sh_on - 1;
    end
    
    sh_off = p_off + 1;
    while sh_off < numel(t) && is_decreasing(sh_off)
        sh_off = sh_off + 1;
    end
    
    if sh_on < 1
        sh_on = 1;
    end
    
    if sh_off > numel(t)
        sh_off = numel(t);
    end
    
    sh(sh_on:sh_off) = true;
end

end
